#Arquivo do proteoma de P.aeruginosa CCBH4851
file= '/home/camila/LMDM/P.Aeruginosa/CCBH4851/CCBH4851-proteome.fasta'
with open(file, 'r') as f:
    #Abre o arquivo no modo de leitura
    for line in f:
        #Percorre as linhas do arquivo
        if line[0]=='>':
            #Se a linha começa com > (um indicador de início de sequência do formato fasta) cria um novo arquivo e escreve nele o conteúdo da linha
            name = line.split()[0].replace('>', '').replace(' ','_') + '.fasta'
            with open(name,'w') as newfile:
                newfile.write(line)
        else:
            #Se a linha não começa com > o arquivo criado anteriormente é aberto e o conteúdo da linha é escrito nele
            with open(name, 'a') as newfile:
                newfile.write(line)
'''
file= open('/home/camila/LMDM/P.Aeruginosa/CCBH4851/teste/CCBH4851-proteome.fasta','r')
newfile= open('a','w')
for line in file.readlines():
    if line[0]=='>':
        newfile.close()
        name = line[1:13] + '.fasta'
        newfile=open(name,'w')
        newfile.write(line)
    else:
            newfile.write(line)
file.close()
'''
