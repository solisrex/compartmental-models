from nicegui import ui
from lark import Lark, Transformer, v_args
from sympy import *


As = [symbols('a'+str(i)) for i in range(6)]
Bs = [symbols('b'+str(i)) for i in range(12)]
Cs = [symbols('c'+str(i)) for i in range(5)]
Ds = [symbols('d'+str(i)) for i in range(5)]

symbols_list = As + Bs + Cs + Ds
symbol_dictionary = dict(zip([str(x) for x in symbols_list], symbols_list))


matrix_1_2 = Matrix([[As[2],Bs[3],Bs[6]],
                    [0,Bs[0]+Bs[7]+Bs[10],-Bs[7]],
                    [0,-Bs[4],Bs[1]+Bs[4]+Bs[11]],
                    [0,Bs[5],Bs[8]]])


matrix_1_3 = Matrix([[As[2],Bs[3],Bs[6],Bs[9]],
                    [0,Bs[0]+Bs[7]+Bs[10],-Bs[7],-Bs[10]],
                    [0,-Bs[4],Bs[1]+Bs[4]+Bs[11],-Bs[11]],
                    [0,-Bs[5],-Bs[8],Bs[2]+Bs[5]+Bs[8]]])

matrix_1_5 = Matrix([[As[2],Bs[3],Bs[6]+Bs[9]],
                    [0,Bs[0]+Bs[7]+Bs[10],-Bs[7]-Bs[10]],
                    [0,-Bs[4],Bs[1]+Bs[4]],
                    [0,-Bs[5],Bs[2]+Bs[5]]])


matrix_2_2 = Matrix([[Bs[3],Bs[6]],
                    [Bs[0]+Bs[7]+Bs[10]+Cs[0],-Bs[7]],
                    [-Bs[4],Bs[1]+Bs[4]+Bs[11]+Cs[0]],
                    [Bs[5],Bs[8]]])


matrix_2_3 = Matrix([[Bs[3],Bs[6],Bs[9]],
                    [Bs[0]+Bs[7]+Bs[10]+Cs[0],-Bs[7],-Bs[10]],
                    [-Bs[4],Bs[1]+Bs[4]+Bs[11]+Cs[0],-Bs[11]],
                    [-Bs[5],-Bs[8],Bs[2]+Bs[5]+Bs[8]+Cs[0]]])



matrix_2_5 = Matrix([[Bs[3],Bs[6]+Bs[9]],
                    [Bs[0]+Bs[7]+Bs[10]+Cs[0],-Bs[10]],
                    [-Bs[4],Bs[1]+Bs[4]+Cs[0]],
                    [-Bs[5],Bs[2]+Bs[5]+Cs[0]]])

matrix_3_2 = Matrix([[As[2],Bs[6],Bs[9]],
                    [As[4],Bs[7],Bs[10]],
                    [0,Bs[1]+Bs[4]+Bs[11],-Bs[11]],
                    [0,-Bs[8],Bs[2]+Bs[5]+Bs[8]]])

matrix_4_2 = Matrix([[Bs[6],Bs[9]],
                    [Bs[7],Bs[10]],
                    [Bs[1]+Bs[4]+Bs[11]+Cs[0],-Bs[11]],
                    [-Bs[8],Bs[2]+Bs[5]+Bs[8]+Cs[0]]])


rational_1 = Matrix([[Bs[0],-Bs[0]-Bs[5],Bs[5]],
                   [-Bs[2]-Bs[4]-Cs[1],Bs[2],Bs[4]],
                   [Bs[1],-As[2]+Bs[3],As[2]-Bs[1]-Bs[3]],
                   [-As[1]+Bs[1],Bs[3],As[1]-Bs[1]-Bs[3]-Cs[1]],
                   [-As[0]+Bs[0]-Bs[2]-Bs[4],As[0]-Bs[0]+Bs[2]-Bs[5]-Cs[1],Bs[4]+Bs[5]]])


rational_2 = Matrix([[Bs[0],-Bs[0]-Bs[5]-Cs[2],Bs[5]],
                   [-Bs[2]-Bs[4]-Cs[1],Bs[2],Bs[4]],
                   [Bs[1],-As[2]+Bs[3],As[2]-Bs[1]-Bs[3]-Cs[2]],
                   [-As[1]+Bs[1],Bs[3],As[1]-Bs[1]-Bs[3]-Cs[1]],
                   [-As[0]+Bs[0]-Bs[2]-Bs[4]-Cs[2],As[0]-Bs[0]+Bs[2]-Bs[5]-Cs[1],Bs[4]+Bs[5]]])



@v_args(inline=True)
class CalculateTree(Transformer):
    
    number = int

    def start(self, *args):
        return list(args)

    def swap_rows(self, i, j):
        row_i = int(i)
        row_j = int(j)
        return ('swap',row_i, row_j)

    def scale_row(self, i, s):
        row_i = int(i)
        scale = s
        return ('scale',row_i, scale)
    
    def add_row(self, i, j):
        row_i = int(i)
        row_j = int(j)
        return ('add',row_i, row_j)

    def sub_row(self, i, j):
        row_i = int(i)
        row_j = int(j)
        return ('sub',row_i, row_j)

    def assign_var(self, x, y):
        var_x = symbol_dictionary[str(x)]
        return ('assign',var_x, y)

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

grammar = open('./matrix.lark').read()

calc_parser = Lark(grammar, parser='lalr',transformer=CalculateTree())
calc = calc_parser.parse


matrix_list = {
    'Matrix 1.2' : matrix_1_2,
    'Matrix 1.3' : matrix_1_3,
    'Matrix 1.5' : matrix_1_5,
    'Matrix 2.2' : matrix_2_2,
    'Matrix 2.4' : matrix_2_3,
    'Matrix 2.5' : matrix_2_5,
    'Matrix 3.2' : matrix_3_2,
    'Matrix 4.2' : matrix_4_2,
    'Rational system 1' : rational_1,
    'Rational system 2' : rational_2
}


selections = list(matrix_list.keys())

current_selection = {'selection' : 'Matrix 1.2'}

def change_matrix(e):
    current_selection.update(selection=e.value)
    left_hand_panel.refresh()
    ui.run_javascript("MathJax.typesetPromise();")

@ui.refreshable
def left_hand_panel():
    transformed_matrix = matrix_list[current_selection['selection']].copy()
    try:
        instructions = calc(form_data['raw_text'])
        if type(instructions) is tuple:
            instructions = [instructions]
        for command in instructions:
            if command[0] == 'swap':
                _, row_i, row_j = command
                transformed_matrix.row_swap(row_i-1, row_j-1)
            elif command[0] == 'scale':
                _, row_i, scale = command
                transformed_matrix.row_op(row_i-1, lambda x, j: scale * x)
            elif command[0] == 'add':
                _, row_i, row_j = command
                transformed_matrix.row_op(row_i-1, lambda x, j: x + transformed_matrix[row_j-1, j])
            elif command[0] == 'sub':
                _, row_i, row_j = command
                transformed_matrix.row_op(row_i-1, lambda x, j: x - transformed_matrix[row_j-1, j])
            elif command[0] == 'assign':
                _, var, value = command
                transformed_matrix = transformed_matrix.subs(var, value)
            transformed_matrix = transformed_matrix.applyfunc(lambda x: simplify(together(x)))
    except Exception as e:
        print(f"Error occurred while parsing: {e}")

    with ui.row().style('width: 100%;'):
        ui.html(r"<b>Original matrix</b><br> $$\Large" + latex(matrix_list[current_selection['selection']],mat_str="matrix") + " $$").style('margin:auto;text-align:center;')
        ui.html(r"<b>Transformed matrix</b><br> $$\Large " + latex(transformed_matrix,mat_str="matrix") + "$$").style('margin:auto;text-align:center;')
        

form_data = {'raw_text': ''}

ui.add_css('''textarea{line-height: 2;font-family: monospace; font-size:20px;height: 80dvh;}''')


text_text = r"\left[\begin{matrix}a_{2} & b_{3} & b_{6}\\ 0 & b_{0} + b_{10} + b_{7} & - b_{7}\\0 & - b_{4} & b_{1} + b_{11} + b_{4}\\0 & b_{5} & b_{8}\end{matrix}\right]"

head_html = '''
<script type="text/javascript" id="MathJax-script" async
  src="https://cdn.jsdelivr.net/npm/mathjax@3/es5/tex-chtml.js">
</script>
<script>
window.addEventListener("load", () => {
    MathJax.typesetPromise();
});
</script>
'''


def update_button():
    left_hand_panel.refresh()
    ui.run_javascript("MathJax.typesetPromise();")

ui.add_head_html(head_html)
with ui.splitter(value=60).style('width: 100%;') as splitter:
    with splitter.before:
        with ui.row().style('margin:auto'):
            ui.label('Matrix Manipulator').style('font-size: 50px; font-weight: bold; margin:auto;')
            ui.select(options=selections,value=selections[0],on_change=lambda e : change_matrix(e))
        left_hand_panel()
        

    with splitter.after:
        ui.textarea('Enter commands:').style('width: 100%;').bind_value(form_data, 'raw_text')
        ui.button('Update',on_click=update_button).style('margin: auto; margin-top: 10px;').classes('col-span-full')

ui.run()
