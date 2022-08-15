import numpy
import warnings

###########################
# Some Import Statements  #
# with exception handling #
###########################

try:
    import tequila as tq
except:
    raise Exception("no tq found, install with pip install tequila-basic")

try:
    from environment import CircuitEnv
except Exception as E:
    print(str(E))
    raise Exception("could not find original code for CircuitEnv, either move this file in the same directory as `` or provide the path over system variable `export PYTHONPATH=$PYTHONPATH:/path/to/directory/with/environment.py`")

##################################################################################
# Derived Class that uses the old encoding to create different types of circuits #
##################################################################################

class CircuitEnvUCC(CircuitEnv):
    """
    Class is derived from CircuitEnv, means it inherits all the functions from CircuitEnv
    The functions which are re-defined here, are overridden (so they will be called instead of theones from CircuitEnv)
    In order to use: in main.py just initialize CircuitEnvUCC instead of CircuitEnv
    E.g. by typing at the top: from ucc import CircuitEnvUCC as CircuitEnv
    with the "as CircuitEnv" part you can avoid to rename the instances in main.py where CircuitEnv is called
    """
    
    def __init__(self, *args, **kwargs):
        """
        Overriding init in order to be able to set the n_electrons variable
        """
        
        # regular initialization
        super().__init__(*args, **kwargs)
        if "n_electrons" in kwargs:
            self.n_electrons = kwargs["n_electrons"]
        else:
            warnings.warn("n_electrons not given ... assuming it's the test case with 2 active electrons", UserWarning)
            self.n_electrons = 2

    def get_energy(self, *args, **kwargs):
        """
        Just illustrating some things on how overriding works
        here we just call the original get_energy function from the base class and return it
        so right now this function dosn't really do anything more and could also be deleted
        For illustration purpose I'll give it an option to remove the energy-shift (wich is not really needed)
        """
        energy = super().get_energy(*args)
        if "remove_shift" in kwargs and kwargs["remove_shift"]:
            energy -= self.energy_shift # energy_shift is set in the initialization of the base class (this is inherited as well)
        return energy

    def make_circuit(self, thetas=None):
        """
        We use the same encoding as in the paper but with different interpretation

        1. We always start with a fixed arrangements of X gates that prepare the Hartree-Fock state
           --> Means this needs to be Jordan-Wigner encoding right now
        2. We use different gates in the encoding (instead or R and CNOT we use single and double qubit excitations: i.e. we are moving "1"s around in the wavefunction)

        The first two integers encode approximations to single-electron rotations in the form of QubitExcitations
        The second two integers encode approximations to double-electron rotatiosn in the form of spin-paired quasi-particles
        Integers encode spatial orbitals which are converted to spin-orbitals (2*x for spin-up, 2*x+1 for spin-down) in Jordan-Wigner encoding)
        i.e. state = [a,b,c,d]

        if c == num_qubits we will initialize a single qubit excitation on spin-up and spin-down electrons
        QubitExcitation(target=[(2*a,2*b)], angle=theta)
        QubitExcitation(target=[(2*a+1,2*b+1)], angle=theta)
        otherwise we will initialize a two-qubit excitation corresponding to an excigtation of a spin-paired pair of electrons
        QubitExcitation(target=[(2*a,2*b,2*a+1,2*b+1)], angle=theta)
        
        Two-Qubit Excitations correspnods to circuit in Eq.(22) of https://arxiv.org/abs/2207.12421
        One-QUbit Excitations corresponds to circuit in Rq.(6) - with more details in the appendix equations (b3),(b4),(b5) of https://arxiv.org/abs/2207.12421 without the Z components (easier for now - can be added later)

        Note for the single excitations: As we are neglecting Z components but at the same time enforcing the same parameters for spin-up and spin-down excitations we are not consistent in the spin symmetry. For LiH this will however not matter much. For later there are two ways to proceed: Either keep the approximation but allow individual angles for spin-up/spin-down excitations in order to variationally counter the neglegt-Z approximation. Or use propper fermionic excitations (i.e. keep the Zs)

        """

        state = self.state.clone()
        U = tq.gates.X([q for q in range(self.n_electrons)])
        for i in range(self.num_layers):
            # at self.reset, make_circuit is called with thetas=state
            # doesn't make much sense so we will initialize only the HF state here
            # see the try/error part below that catches this 
            a = int(state[0][i].item()//2) 
            b = int(state[1][i].item()//2)
            c = int(state[2][i].item()//2)
            d = int(state[3][i].item()//2)
            try:
                angle = thetas[i].item()*numpy.pi
                # angle needs to be set
                assert angle is not None
            except:
                continue

            if c != self.num_qubits//2 and c!=d:
                U += tq.gates.QubitExcitation(target=[2*c,2*d], angle=angle)
                U += tq.gates.QubitExcitation(target=[2*c+1,2*d+1], angle=angle)
            elif a != self.num_qubits//2 and a!=b:
                U += tq.gates.QubitExcitation(target=[2*a,2*b,2*a+1,2*b+1], angle=angle)
            
        
        # make sure that tq doesn't do any automatic mappings to smaller qubit systems
        U.n_qubits=self.num_qubits
        U += tq.gates.X([i for i in range(self.num_qubits)])
        U += tq.gates.X([i for i in range(self.num_qubits)])
        # failsave as the RL code works with explicit variables
        assert len(U.extract_variables()) == 0 
        # convert to Qulacs circuit
        circuit = tq.compile(U, backend="qulacs").circuit
    
        return circuit
