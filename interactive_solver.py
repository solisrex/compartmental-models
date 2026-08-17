from nicegui import ui
from lark import Lark, Transformer, v_args
from sympy import Not, Or, re, im, I, E, symbols, latex, Eq, Function, simplify

As = [symbols('a'+str(i)) for i in range(6)]
Bs = [symbols('b'+str(i)) for i in range(12)]
Cs = [symbols('c'+str(i)) for i in range(5)]
Ds = [symbols('d'+str(i)) for i in range(5)]


X1 = Function('X1')(symbols('t'))
X2 = Function('X2')(symbols('t'))   
X3 = Function('X3')(symbols('t'))
X4 = Function('X4')(symbols('t'))

X1_dot =  As[0]*X1*X2 + As[1]*X1*X3 + As[2]*X1*X4 - (Bs[3]+Bs[6]+Bs[9])*X1 + Bs[0]*X2 + Bs[1]*X3 + Bs[2]*X4;

X2_dot = -As[0]*X1*X2 + As[3]*X2*X3 + As[4]*X2*X4 + Bs[3]*X1 - (Bs[0]+Bs[7]+Bs[10])*X2 + Bs[4]*X3 + Bs[5]*X4;

X3_dot = -As[1]*X1*X3 - As[3]*X2*X3 + As[5]*X3*X4 + Bs[6]*X1 + Bs[7]*X2 - (Bs[1]+Bs[4]+Bs[11])*X3 + Bs[8]*X4;

X4_dot = -As[2]*X1*X4 - As[4]*X2*X4 - As[5]*X3*X4 + Bs[9]*X1 + Bs[10]*X2 + Bs[11]*X3 - (Bs[2]+Bs[5]+Bs[8])*X4;

darboux = Ds[0] + Ds[1]*X1 + Ds[2]*X2 + Ds[3]*X3 + Ds[4]*X4

cofactor = Cs[0] + Cs[1]*X1 + Cs[2]*X2 + Cs[3]*X3 + Cs[4]*X4

derivatives = [X1_dot, X2_dot, X3_dot, X4_dot]

derivative_labels = ['X_{1}\'', 'X_{2}\'', 'X_{3}\'', 'X_{4}\'']

darboux_label = 'D'
cofactor_label = 'C'

eqns = [
    Eq(Cs[0]*Ds[0], 0),
    Eq(Cs[1]*Ds[1], 0),
    Eq(Cs[2]*Ds[2], 0),
    Eq(Cs[3]*Ds[3], 0),
    Eq(Cs[4]*Ds[4], 0),
    Eq(As[0]*Ds[1] - Cs[2]*Ds[1] - As[0]*Ds[2] - Cs[1]*Ds[2], 0),
    Eq(As[1]*Ds[1] - Cs[3]*Ds[1] - As[1]*Ds[3] - Cs[1]*Ds[3], 0),
    Eq(As[2]*Ds[1] - Cs[4]*Ds[1] - As[2]*Ds[4] - Cs[1]*Ds[4], 0),
    Eq(As[3]*Ds[2] - Cs[3]*Ds[2] - As[3]*Ds[3] - Cs[2]*Ds[3], 0),
    Eq(As[4]*Ds[2] - Cs[4]*Ds[2] - As[4]*Ds[4] - Cs[2]*Ds[4], 0),
    Eq(As[5]*Ds[3] - Cs[4]*Ds[3] - As[5]*Ds[4] - Cs[3]*Ds[4], 0),
    Eq(-Cs[1]*Ds[0] - Bs[3]*Ds[1] - Bs[6]*Ds[1] - Bs[9]*Ds[1] - Cs[0]*Ds[1] + Bs[3]*Ds[2] + Bs[6]*Ds[3] + Bs[9]*Ds[4], 0),
    Eq(-Cs[2]*Ds[0] + Bs[0]*Ds[1] - Bs[0]*Ds[2] - Bs[10]*Ds[2] - Bs[7]*Ds[2] - Cs[0]*Ds[2] + Bs[7]*Ds[3] + Bs[10]*Ds[4], 0),
    Eq(-Cs[3]*Ds[0] + Bs[1]*Ds[1] + Bs[4]*Ds[2] - Bs[1]*Ds[3] - Bs[11]*Ds[3] - Bs[4]*Ds[3] - Cs[0]*Ds[3] + Bs[11]*Ds[4], 0),
    Eq(-Cs[4]*Ds[0] + Bs[2]*Ds[1] + Bs[5]*Ds[2] + Bs[8]*Ds[3] - Bs[2]*Ds[4] - Bs[5]*Ds[4] - Bs[8]*Ds[4] - Cs[0]*Ds[4], 0)
]

symbols_list = As + Bs + Cs + Ds

symbol_dictionary = dict(zip([str(x) for x in symbols_list], symbols_list))

@v_args(inline=True)
class CalculateTree(Transformer):
    
    number = int

    def start(self, *args):
        return list(args)

    def assign_var(self, x, y):
        var_x = symbol_dictionary[str(x)]
        return (var_x, y)

    def var(self, name):
        return symbol_dictionary[str(name)]

    def __init__(self):
        self.vars = {}

    def neg(self, x):
        return -x

    def add(self, x, y):
        return x + y
    
    def sub(self, x, y):
        return x - y
    
    def mul(self, x, y):
        return x * y
    
    def div(self, x, y):
        return x / y

    def pow(self, x, y):
        return x ** y

    def test(self, *args):
        print('test called with:', args)

grammar = open('solver.lark').read()

calc_parser = Lark(grammar, parser='lalr',transformer=CalculateTree())
calc = calc_parser.parse


system_text = '''
$$
X_{1}' =  a_{0} X_{1} X_{2} + a_{1} X_{1} X_{3} + a_{2} X_{1} X_{4} - (b_{3}+b_{6}+b_{9}) X_{1} + b_{0} X_{2} + b_{1} X_{3} + b_{2} X_{4};\\\\
$$

$$
X_{2}' = -a_{0} X_{1} X_{2} + a_{2} X_{2} X_{3} + a_{3} X_{2} X_{4} + b_{3} X_{1} - (b_{0}+b_{7}+b_{10}) X_{2} + b_{4} X_{3} + b_{5} X_{4};\\\\
$$

$$
X_{3}' = -a_{1} X_{1} X_{3} - a_{2} X_{2} X_{3} + a_{4} X_{3} X_{4} + b_{6} X_{1} + b_{7} X_{2} - (b_{1}+b_{4}+b_{11}) X_{3} + b_{8} X_{4};\\\\
$$

$$
X_{4}' = -a_{2} X_{1} X_{4} - a_{3} X_{2} X_{4} - a_{4} X_{3} X_{4} + b_{9} X_{1} + b_{10} X_{2} + b_{11} X_{3} - (b_{2}+b_{5}+b_{8}) X_{4};
$$
'''

@ui.refreshable
def left_hand_panel():
    try:
        substitutions = calc(form_data['raw_text'])
        if type(substitutions) is not list:
            substitutions = [substitutions]
    except:
        print(f"Error occurred while parsing.")
        substitutions = []

    ui.label('System of Equations').style('margin: auto;font-size: 20px;').classes('col-span-full')
    for eqn,label in zip(derivatives + [darboux,cofactor], derivative_labels + [darboux_label, cofactor_label]):
        subbed_item = eqn.subs(substitutions)
        ui.markdown('$$'+label+'='+latex(subbed_item)+'$$', extras=['latex']).style('font-size: 15px;')
    ui.label('Darboux polynomial').style('margin: auto;font-size: 20px;').classes('col-span-full')
    
    with ui.scroll_area().style('height: 60dvh;'):
        with ui.grid(columns='1fr',).style('align-items: center; justify-content:right;width: 100%;'):
            for eqn in eqns:
                subbed_eqn = simplify(eqn.subs(substitutions))
                
                if subbed_eqn == True:
                    ui.label('✅')
                elif subbed_eqn == False:
                    ui.label('❌')
                else:
                    ui.markdown('$'+latex(subbed_eqn)+'$', extras=['latex']).style('font-size: 20px;')

                ui.separator().classes('col-span-full')


form_data = {'raw_text': ''}

ui.add_css('''textarea{line-height: 2;font-family: monospace; font-size:20px;height: 80dvh;}''')

with ui.splitter(value=60).style('width: 100%;') as splitter:
    with splitter.before:
        left_hand_panel()


    with splitter.after:
        
        ui.label('Variable values').style('margin: auto;font-size: 20px;').classes('col-span-full')
        ui.textarea('Assign values to variables:').style('width: 100%;').bind_value(form_data, 'raw_text')
        ui.button('Update',on_click=lambda x : left_hand_panel.refresh()).style('margin: auto; margin-top: 10px;').classes('col-span-full')

ui.run()
