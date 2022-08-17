import numpy
import warnings
import torch
import copy

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

        if "geometry" in kwargs:
            geometry = kwargs["geometry"]
        else:
            geometry = "Li 0.0 0.0 0.0\nH 0.0 0.0 2.2"

        if "active_orbitals" in kwargs:
            active_orbitals = kwargs["active_orbitals"]
        elif "li" in geometry.lower() and "h" in geometry.lower():
            active_orbitals = [1, 2, 5]
        else:
            active_orbitals = None

        self.mol = tq.Molecule(geometry=geometry, basis_set="sto-3g", active_orbitals=active_orbitals)
        self.tq_hamiltonian = self.mol.make_hamiltonian()
        assert self.n_electrons == self.mol.n_electrons

        # consistency check

        __ham = numpy.load(f"mol_data/LiH_{self.num_qubits}q_geom_{self.geometry}_{self.ham_mapping}.npz")
        self.hamiltonian, eigvals, self.energy_shift = __ham['hamiltonian'], __ham['eigvals'], __ham['energy_shift']

        min_eig = min(eigvals)+self.energy_shift
        max_eig = max(eigvals)+self.energy_shift
        if self.tq_hamiltonian.n_qubits < 11:
            v = numpy.linalg.eigvalsh(self.tq_hamiltonian.to_matrix())

            consistent = True
            if not numpy.isclose(v[-1], max_eig):
                consistent = False
            if not numpy.isclose(v[0],min_eig):
                consistent = False

            if not consistent:
                warnings.warn("tq molecule and loaded qiskit Hamiltonian have different energies!")

        if self.num_layers > 15:
            warnings.warn("num_layers is quite high", UserWarning)



    def get_energy(self, *args, **kwargs):
        """
        Make sure this is not called
        """
        raise Exception("forbidden")

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
        variables = {}
        U = tq.gates.X([q for q in range(self.n_electrons)])
        for i in range(self.num_layers):
            # at self.reset, make_circuit is called with thetas=state
            # doesn't make much sense so we will initialize only the HF state here
            # see the try/error part below that catches this 
            a = int(state[0][i].item()//2) 
            b = int(state[1][i].item()//2)
            c = int(state[2][i].item()//2)
            d = int(state[3][i].item()//2)

            angle = None
            if a < self.num_qubits//2 and b < self.num_qubits//2 and a!=b:
                angle=(a,b,"D")
                U += tq.gates.QubitExcitation(target=[2*a,2*b,2*a+1,2*b+1], angle=tq.Variable(angle)*numpy.pi)
            elif c != self.num_qubits//2 and d < self.num_qubits//2 and c!=d:
                angle=(c,d,"S")
                U += tq.gates.QubitExcitation(target=[2*c,2*d], angle=(tq.Variable(angle)+0.2)*numpy.pi)
                U += tq.gates.QubitExcitation(target=[2*c+1,2*d+1], angle=(tq.Variable(angle)+0.2)*numpy.pi)

            if angle is not None:
                try:
                    variables[angle] = thetas[i].item()
                except:
                    variables[angle] = 0.0

        self.variables = variables
        return U

    def step(self, action, train_flag=True):

        """
        Action is performed on the first empty layer.
        Variable 'actual_layer' points last non-empty layer.
        """
        next_state = self.state.clone()
        self.actual_layer += 1

        """
        First two elements of the 'action' vector describes position of the CNOT gate.
        Position of rotation gate and its axis are described by action[2] and action[3].
        When action[0] == num_qubits, then there is no CNOT gate.
        When action[2] == num_qubits, then there is no Rotation gate.
        """

        next_state[0][self.actual_layer] = action[0]
        next_state[1][self.actual_layer] = (action[0] + action[1]) % self.num_qubits

        ## state[2] corresponds to number of qubit for rotation gate
        next_state[2][self.actual_layer] = action[2]
        next_state[3][self.actual_layer] = action[3]
        next_state[4][self.actual_layer] = torch.zeros(1)

        self.state = next_state.clone()

        # solve with scipyt and analytical gradients (will take longer)
        U = self.make_circuit()
        H = self.tq_hamiltonian
        E = tq.ExpectationValue(H=H, U=U)
        result = tq.minimize(E, silent=True, initial_values=self.variables)
        energy = result.energy
        angles = list(result.variables.values())
        thetas = self.state[-1]
        for i in range(len(angles)):
            thetas[i] = angles[i]
        next_state[-1] = thetas

        self.energy = energy
        if energy < self.curriculum.lowest_energy and train_flag:
            self.curriculum.lowest_energy = copy.copy(energy)

        self.error = float(abs(self.min_eig - energy))

        rwd = self.reward_fn(energy)
        self.prev_energy = numpy.copy(energy)

        energy_done = int(self.error < self.done_threshold)
        layers_done = self.actual_layer == (self.num_layers - 1)
        done = int(energy_done or layers_done)

        if done:
            self.curriculum.update_threshold(energy_done=energy_done)
            self.done_threshold = self.curriculum.get_current_threshold()
            self.curriculum_dict[str(self.current_bond_distance)] = copy.deepcopy(self.curriculum)

        self.state = next_state.clone()
        return next_state.view(-1).to(self.device), torch.tensor(rwd, dtype=torch.float32, device=self.device), done


