# dependencies:
# tequila: pip install tequila-basic
# qpic: pip install qpic
# working latex compiler
# working pdf to png command
# if the latter is not present: replace png by pdf
import tequila as tq

# pic for single qubit excitation
U = tq.gates.QubitExcitation(target=[0,2],angle="a")
U+= tq.gates.QubitExcitation(target=[1,3],angle="a")
U = tq.compile_circuit(U)
U.export_to(filename="single_excitation.png")

# pic for double qubit excitation
U = tq.gates.QubitExcitation(target=[0,2,1,3],angle="a")
U = tq.compile_circuit(U)
U.export_to(filename="double_excitation.png")

