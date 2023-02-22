#  X = TAMANHO(abcissa)
#  Y = VALOR TESTE(ordenada)
import matplotlib.pyplot as plt
import pandas as pd

def grep(file_path, substring):
    myLine = {'line': '', 'index':0}
    with open(file_path, 'r') as file:
        i = 0
        for line in file:
            if substring in line:
                myLine['line'] =line 
                myLine['index'] =i  
                i+=1
                yield myLine

TAMANHOS=['32','64', '128', '256', '512', '1000', '2000', '4000', '8000']
HEADER = ['OP1 v1','OP2 v1','OP1 v2','OP2 v2']

TYPESDICTIONARY = {
    'L3': 'L3 bandwidth [MBytes/s]', 
    'L2CACHE': 'L2 miss ratio'
}
    # 'time':'RDTSC Runtime [s]'
    # 'FLOPS_DP' : 'DP MFLOP/s'
CSVPATH = 'CSV/'
IMGSPATH = 'IMGS_PLOTS/'

def generate_csv(file_path, type, typeMetric, typeHeader):
    versoes = ['v1', 'v2']
    operacoes =['op1','op2']

    with open(file_path, 'w') as f:
                print('Tamanho', end=" ",file=f)
                for col in typeHeader:
                    print(",",col, end=" ",file=f)
                print(file=f)
                for tam in TAMANHOS:
                    print(tam, end=" ",file=f)
                    for myLine in grep('./v1'+'/logs/lik'+tam+type+'.log', typeMetric):
                        formattedLine = myLine['line'].replace(',', ' ') 
                        words = formattedLine.split()
                        print(",",words[-1], end=" ", file=f)
                
                    for myLine in grep('./v2'+'/logs/lik'+tam+type+'.log', typeMetric):
                        formattedLine = myLine['line'].replace(',', ' ')
                        words = formattedLine.split()
                        print(",",words[-1], end=" ", file=f)
                    
                    print(file=f)

def plot_multiline_result(input_file_path, title, output_file_path, ytitle):
    df = pd.read_csv(input_file_path)

    for i in range(1, len(df.columns), 2):
        plt.plot(df[df.columns[0]], df[df.columns[i]],marker='v',label=df.columns[i], linestyle='-',markersize=2, linewidth=1) # colunas nao otimizadas
        plt.plot(df[df.columns[0]], df[df.columns[i+1]],marker='o',label=df.columns[i+1], linestyle='--',markersize=2, linewidth=1) # colunas otimizadas

    plt.xlabel(df.columns[0])
    plt.ylabel(ytitle)
    plt.title(title)
    plt.legend(df.columns[1:])

    plt.minorticks_on()
    plt.grid(which='minor', axis='y', linestyle=':', linewidth=0.5)
    # plt.show()
    plt.savefig(output_file_path, dpi=300)
    plt.clf()

for type in TYPESDICTIONARY.keys():
    # generate_csv(CSVPATH+type+'.csv', type, TYPESDICTIONARY[type], HEADER)
    plot_multiline_result(CSVPATH+type+'.csv', type, IMGSPATH+type+'.png',TYPESDICTIONARY[type])

# generate_csv(CSVPATH+'FLOPS_DP.csv', 'FLOPS_DP', 'DP MFLOP/s', HEADER)
plot_multiline_result(CSVPATH+'FLOPS_DP.csv', 'FLOPS_DP',IMGSPATH+ 'FLOPS_DP.png','DP MFLOP/s')
# generate_csv(CSVPATH+'TIME.csv', 'L3', 'Runtime (RDTSC)', HEADER)
plot_multiline_result(CSVPATH+'TIME.csv', 'TIME',IMGSPATH+ 'time.png','Runtime (RDTSC) s')

# Banda de Memória: MEM ou L3₢C
# Memory bandwidth [MBytes/s] 
# Cache miss L2: CACHE ou L2CACHE
# data cache miss ratio
# Operações aritméticas: FLOPS_DP 
# FLOPS_DP  e FLOPS_AVX 