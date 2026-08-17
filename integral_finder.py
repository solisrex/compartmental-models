from sympy import solve, symbols, Function, Eq, latex

As = [symbols('a'+str(i)) for i in range(6)]
Bs = [symbols('b'+str(i)) for i in range(12)]

X1 = Function('X1')(symbols('t'))
X2 = Function('X2')(symbols('t'))   
X3 = Function('X3')(symbols('t'))
X4 = Function('X4')(symbols('t'))

X1_dot =  As[0]*X1*X2 + As[1]*X1*X3 + As[2]*X1*X4 - (Bs[3]+Bs[6]+Bs[9])*X1 + Bs[0]*X2 + Bs[1]*X3 + Bs[2]*X4;
X2_dot = -As[0]*X1*X2 + As[3]*X2*X3 + As[4]*X2*X4 + Bs[3]*X1 - (Bs[0]+Bs[7]+Bs[10])*X2 + Bs[4]*X3 + Bs[5]*X4;
X3_dot = -As[1]*X1*X3 - As[3]*X2*X3 + As[5]*X3*X4 + Bs[6]*X1 + Bs[7]*X2 - (Bs[1]+Bs[4]+Bs[11])*X3 + Bs[8]*X4;
X4_dot = -As[2]*X1*X4 - As[4]*X2*X4 - As[5]*X3*X4 + Bs[9]*X1 + Bs[10]*X2 + Bs[11]*X3 - (Bs[2]+Bs[5]+Bs[8])*X4;

derivatives = [X1_dot,X2_dot,X3_dot,X4_dot]
derivative_labels = ['\\dot{X}_1','\\dot{X}_2','\\dot{X}_3','\\dot{X}_4']


trans12 = {
    X1:X2,
    X2:X1,
    As[0]:-As[0],
    As[1]:As[3],
    As[3]:As[1],
    As[2]:As[4],
    As[4]:As[2],
    Bs[0]:Bs[3],
    Bs[3]:Bs[0],
    Bs[6]:Bs[7],
    Bs[7]:Bs[6],
    Bs[9]:Bs[10],
    Bs[10]:Bs[9],
    Bs[1]:Bs[4],
    Bs[4]:Bs[1],
    Bs[2]:Bs[5],
    Bs[5]:Bs[2]}

trans13 = {
    X1:X3,
    X3:X1,
    As[1]:-As[1],
    As[0]:-As[3],
    As[3]:-As[0],
    As[2]:As[5],
    As[5]:As[2],
    Bs[0]:Bs[7],
    Bs[7]:Bs[0],
    Bs[6]:Bs[1],
    Bs[1]:Bs[6],
    Bs[2]:Bs[8],
    Bs[8]:Bs[2],
    Bs[3]:Bs[4],
    Bs[4]:Bs[3],
    Bs[9]:Bs[11],
    Bs[11]:Bs[9]}

trans14 = {
    X1:X4,
    X4:X1,
    As[2]:-As[2],
    As[0]:-As[4],
    As[4]:-As[0],
    As[1]:-As[5],
    As[5]:-As[1],
    Bs[2]:Bs[9],
    Bs[9]:Bs[2],
    Bs[0]:Bs[10],
    Bs[10]:Bs[0],
    Bs[1]:Bs[11],
    Bs[11]:Bs[1],
    Bs[3]:Bs[5],
    Bs[5]:Bs[3],
    Bs[8]:Bs[6],
    Bs[6]:Bs[8]}

trans23 = {
    X2:X3,
    X3:X2,
    As[0]:As[1],
    As[1]:As[0],
    As[3]:-As[3],
    As[4]:As[5],
    As[5]:As[4],
    Bs[3]:Bs[6],
    Bs[6]:Bs[3],
    Bs[4]:Bs[7],
    Bs[7]:Bs[4],
    Bs[5]:Bs[8],
    Bs[8]:Bs[5],
    Bs[0]:Bs[1],
    Bs[1]:Bs[0],
    Bs[10]:Bs[11],
    Bs[11]:Bs[10]}

trans24 = {
    X2:X4,
    X4:X2,
    As[0]:As[2],
    As[2]:As[0],
    As[4]:-As[4],
    As[3]:-As[5],
    As[5]:-As[3],
    Bs[3]:Bs[9],
    Bs[9]:Bs[3],
    Bs[5]:Bs[10],
    Bs[10]:Bs[5],
    Bs[4]:Bs[11],
    Bs[11]:Bs[4],
    Bs[0]:Bs[2],
    Bs[2]:Bs[0],
    Bs[7]:Bs[8],
    Bs[8]:Bs[7]}

trans34 = {
    X3:X4,
    X4:X3,
    As[1]:As[2],
    As[2]:As[1],
    As[3]:As[4],
    As[4]:As[3],
    As[5]:-As[5],
    Bs[6]:Bs[9],
    Bs[9]:Bs[6],
    Bs[7]:Bs[10],
    Bs[10]:Bs[7],
    Bs[8]:Bs[11],
    Bs[11]:Bs[8],
    Bs[2]:Bs[1],
    Bs[1]:Bs[2],
    Bs[4]:Bs[5],
    Bs[5]:Bs[4]}


class Permutation:
    def __init__(self, name, substitution_list):
        self.name = name
        self.substitution_list = substitution_list

    def apply(self, expr):
        permuted_expr = expr
        for substitution in self.substitution_list[::-1]:
            permuted_expr = permuted_expr.subs(substitution,simultaneous=True)
        return permuted_expr


identity = Permutation(name='id', substitution_list=[])
s12 = Permutation(name='(12)', substitution_list=[trans12])
s13 = Permutation(name='(13)', substitution_list=[trans13])
s14 = Permutation(name='(14)', substitution_list=[trans14])
s23 = Permutation(name='(23)', substitution_list=[trans23])
s24 = Permutation(name='(24)', substitution_list=[trans24])
s34 = Permutation(name='(34)', substitution_list=[trans34])
s12_34 = Permutation(name='(12)(34)', substitution_list=[trans34,trans12])
s13_24 = Permutation(name='(13)(24)', substitution_list=[trans24,trans13])
s14_23 = Permutation(name='(14)(23)', substitution_list=[trans23,trans14])
s123 = Permutation(name='(123)', substitution_list=[trans13,trans12]) # (a b c) = (a c)(a b)
s132 = Permutation(name='(132)', substitution_list=[trans12,trans13])
s124 = Permutation(name='(124)', substitution_list=[trans14,trans12])
s142 = Permutation(name='(142)', substitution_list=[trans12,trans14])
s134 = Permutation(name='(134)', substitution_list=[trans14,trans13])
s143 = Permutation(name='(143)', substitution_list=[trans13,trans14])
s234 = Permutation(name='(234)', substitution_list=[trans24,trans23])
s243 = Permutation(name='(243)', substitution_list=[trans23,trans24])
s1234 = Permutation(name='(1234)', substitution_list=[trans14,trans13,trans12]) # (a b c d) = (a d)(a c)(a b)
s1243 = Permutation(name='(1243)', substitution_list=[trans13,trans14,trans12])
s1324 = Permutation(name='(1324)', substitution_list=[trans14,trans12,trans13])
s1342 = Permutation(name='(1342)', substitution_list=[trans12,trans14,trans13])
s1423 = Permutation(name='(1423)', substitution_list=[trans13,trans12,trans14])
s1432 = Permutation(name='(1432)', substitution_list=[trans12,trans13,trans14])

permutations = [identity,s12,s13,s14,s23,s24,s34,s12_34,s13_24,s14_23,
                s123,s132,s124,s142,s134,s143,s234,s243,
                s1234,s1243,s1324,s1342,s1423,s1432]


permutation_names = {i:per.name for i, per in enumerate(permutations)}


class EquivalenceClass:
    def __init__(self,name,constraints,darboux,cofactor):
        self.name = name
        self.constraints = constraints
        self.darboux = darboux
        self.cofactor = cofactor

equiv_class_1_1 = EquivalenceClass(
    name='Equivalence class 1.1',
    constraints = [Eq(As[3],0),Eq(As[4],0),Eq(Bs[4],0),Eq(Bs[5],0),Eq(Bs[0]+Bs[7]+Bs[10],0)],
    darboux = Bs[3] - As[0]*X2,
    cofactor = -As[2]*X1)

equiv_class_1_4 = EquivalenceClass(
    name='Equivalence class 1.4',
    constraints = [Eq(As[0],As[1]),Eq(As[4],0),Eq(As[5],0),
                   Eq(Bs[0]+Bs[10],0),Eq(Bs[1]+Bs[11],0),Eq(Bs[5]+Bs[8],0)],
    darboux = (Bs[3]+Bs[6]) - As[1]*(X2+X3),
    cofactor = -As[2]*X1)


equiv_class_1_6 = EquivalenceClass(
    name='Equivalence class 1.6',
    constraints = [Eq(As[0],As[1]),Eq(As[1],As[2]),Eq(Bs[0],0),Eq(Bs[1],0),Eq(Bs[2],0)],
    darboux = (Bs[3]+Bs[6]+Bs[9]) - As[2]*(X2+X3+X4),
    cofactor = -As[0]*X1)

equiv_class_2_1 = EquivalenceClass(
    name='Equivalence class 2.1',
    constraints = [Eq(As[3],0),Eq(As[4],0),Eq(Bs[3],0),Eq(Bs[4],0),Eq(Bs[5],0)],
    darboux = X2,
    cofactor = -(Bs[0]+Bs[7]+Bs[10])-As[0]*X1)

equiv_class_2_4 = EquivalenceClass(
    name='Equivalence class 2.4',
    constraints = [Eq(As[0],As[1]),Eq(As[4],0),Eq(As[5],0),
                   Eq(Bs[4]+Bs[6],0),Eq(Bs[5]+Bs[8],0),
                   Eq(Bs[0]+Bs[10]-Bs[1]-Bs[11],0)],
    darboux = X2+X3,
    cofactor = -(Bs[0]+Bs[10])-As[1]*X1)

equiv_class_2_6 = EquivalenceClass(
    name='Equivalence class 2.6',
    constraints = [Eq(As[0],As[1]),Eq(As[1],As[2]),
                   Eq(Bs[0],Bs[1]),Eq(Bs[1],Bs[2]),
                   Eq(As[3],0),Eq(Bs[3]+Bs[6]+Bs[9],0)],
    darboux = X2+X3+X4,
    cofactor = -Bs[0]-As[2]*X1)

equiv_class_4_1 = EquivalenceClass(
    name='Equivalence class 4.1',
    constraints = [Eq(As[5],0),Eq(Bs[6],0),
                   Eq(Bs[7],0),Eq(Bs[8],0)],
    darboux = X3,
    cofactor = -(Bs[1]+Bs[4]+Bs[11])-As[1]*X1-As[3]*X2)

equiv_class_4_3 = EquivalenceClass(
    name='Equivalence class 4.3',
    constraints = [Eq(As[1],As[2]),Eq(As[3],As[4]),
                   Eq(Bs[6]+Bs[9],Bs[7]+Bs[10])],
    darboux = X3+X4,
    cofactor = -(Bs[6]+Bs[9])-As[2]*X1-As[4]*X2)

equivalence_classes = [equiv_class_1_1,equiv_class_1_4,equiv_class_1_6,
                       equiv_class_2_1,equiv_class_2_4,equiv_class_2_6,
                       equiv_class_4_1,equiv_class_4_3]

equvialence_class_names = {i:equiv_class.name for i, equiv_class in enumerate(equivalence_classes)}


from nicegui import ui

darboux_cards = [(0,0),(1,0)]

def darboux_card(darboux_index,class_index,permutation_index):
    @ui.refreshable
    def display_polynomials(class_index,permutation_index):
        with ui.row():
            permutation = permutations[permutation_index]
            equiv_class = equivalence_classes[class_index]
            permuted_constrains = [permutation.apply(constraint) for constraint in equiv_class.constraints]
            constraint_string = ', '.join(['$'+latex(constraint)+'$' for constraint in permuted_constrains]) 
            permuted_darboux = permutation.apply(equiv_class.darboux)
            permuted_cofactor = permutation.apply(equiv_class.cofactor)                   
            ui.markdown(f"Constraints: {constraint_string}<br>Darboux: ${latex(permuted_darboux)}$<br>Cofactor: ${latex(permuted_cofactor)}$",extras=['latex'])

    def update_class_index(e):
        darboux_cards[darboux_index] = (e.value,permutation_index)
        display_cards.refresh()

    def update_permutation_index(e):
        darboux_cards[darboux_index] = (class_index,e.value)
        display_cards.refresh()

    def delete_darboux(e):
        darboux_cards.pop(darboux_index)
        display_cards.refresh()

    with ui.card().style('margin:auto;width:90%;'):
        ui.html(f'<h4>Darboux {darboux_index + 1}</h4>')
        ui.button(icon='close',on_click=delete_darboux).props('fab color=accent').style("position:absolute; top:0;right: 0;transform: translate(100%,-100%);scale:50%")
        with ui.row():
            ui.label("Equivalence class").style('margin:auto')
            select_class = ui.select(equvialence_class_names,value=class_index,on_change=update_class_index).style('margin:auto')
            ui.label("Permutation").style('margin:auto')
            select_permutation = ui.select(permutation_names,value=permutation_index,on_change=update_permutation_index).style('margin:auto')
        ui.separator()
        display_polynomials(class_index,permutation_index)

@ui.refreshable
def display_cards():
    if not darboux_cards:
        ui.label('No darbouxs')
        return
    for i, (class_index, permutation_index) in enumerate(darboux_cards):
        darboux_card(i, class_index, permutation_index)

def add_darboux():
    darboux_cards.append((0,0))
    display_cards.refresh()

def merge_darbouxs():
    total_constraints = []
    for class_index, permutation_index in darboux_cards:
        permutation = permutations[permutation_index]
        equiv_class = equivalence_classes[class_index]
        permuted_constrains = [permutation.apply(constraint) for constraint in equiv_class.constraints]
        total_constraints.extend(permuted_constrains)
    solutions = solve(total_constraints,As+Bs)
    system_of_equations.refresh(solutions)
    ui.run_javascript("MathJax.typesetPromise();")


@ui.refreshable
def system_of_equations(solutions):

    with ui.row():
        eqn_string = r'\\'.join([f'{eqn_label}&={latex(eqn.subs(solutions))}' for eqn_label, eqn in zip(derivative_labels,derivatives)])
        eqn_start = r'\begin{align*}'
        eqn_end = r'\end{align*}'
        ui.html(f"<h3>System of equations</h3>{eqn_start}{eqn_string}{eqn_end}").style('scale:0.9;')
   
    with ui.row():
        solution_string = r',\,'.join([latex(lhs)+'='+latex(rhs) for lhs,rhs in solutions.items()])
        ui.html(f"<h5>Combined constraints</h5> $${solution_string}$$").style('width:100%')

        ui.html(f"<h5>Darboux/cofactor pairs</h5>").style('width:100%')
        for i, darboux_card in enumerate(darboux_cards):
            class_index, permutation_index = darboux_card
            permutation = permutations[permutation_index]
            equiv_class = equivalence_classes[class_index]
            permuted_darboux = permutation.apply(equiv_class.darboux)
            permuted_cofactor = permutation.apply(equiv_class.cofactor)                   

            ui.html(f"$$D_{{{i+1}}}={latex(permuted_darboux.subs(solutions))},\\quad C_{{{i+1}}}={latex(permuted_cofactor.subs(solutions))}$$").style('width:100%')
                    


total_constraints = []
for class_index, permutation_index in darboux_cards:
    permutation = permutations[permutation_index]
    equiv_class = equivalence_classes[class_index]
    permuted_constrains = [permutation.apply(constraint) for constraint in equiv_class.constraints]
    total_constraints.extend(permuted_constrains)

solutions = solve(total_constraints,As+Bs)




def root():
    head_html = '''
<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js">
</script>
'''
    ui.add_head_html(head_html)
    ui.run_javascript("MathJax.typesetPromise();")
    with ui.row().classes('w-full'):
        with ui.column().classes('flex-1'):
            ui.html('<h3>Integral finder</h3>')
            with ui.column().style('padding-top:50px;padding-bottom:50px;max-height:65dvh; width:100%;overflow-x: hidden;overflow-y: auto;'):
                display_cards()
                ui.button(icon='add',on_click=add_darboux).props('fab color=accent').style("margin:auto; margin-top: 10px;")

            with ui.card().style('width:90%;margin:auto;'):
                ui.textarea(placeholder='Enter additional constraints...').props('autogrow').style('width:100%;')
            ui.button('Merge',on_click=merge_darbouxs).props('fab color=accent').style("margin:auto")

        with ui.column().classes('flex-1'):
            system_of_equations(solutions)
ui.run(root)
