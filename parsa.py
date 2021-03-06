import pandas as pd


with open('dat/parsa.txt') as f:
    line  = f.readline()
    line  = f.readline()
    line  = line.replace('# ', '')
    
    names = line.split()
    
data = pd.read_csv('dat/parsa.txt', sep='\s+', comment ='#', names=names)
data = data.to_latex(index=False)

with open('dat/parsa.tex', 'w') as tf:
     tf.write(data)
