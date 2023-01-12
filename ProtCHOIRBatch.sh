#!/bin/bash
#Estabelece como diretório base o diretório em que o arquivo se encontra
basedir=`pwd`
#Percorre todos os arquivos no formato fasta do diretório base
for file in `ls *.fasta`; do
  oln=`echo $file | awk -F'.' '{print $1}'`
  #Se já existir um arquivo com o final _ProtCHOIR_OUT.zip dentro de uma pasta com o nome do arquivo fasta passa para o próxima interação pois a modelagem para aquele arquivo já foi feita
  if [ -f "${oln}/${oln}_ProtCHOIR_OUT.zip" ]; then
    echo "${oln}/${oln}_ProtCHOIR_OUT.zip exists. Skipping."
#Caso contrário, cria um diretório com o nome do arquivo fasta e executa o programa ProtCHOIR utilizando como entrada o arquivo fasta dentro do novo diretório criado
  else
    mkdir $oln
    cp ${file} ${oln}/${oln}.fasta
    cd ${oln}
    ProtCHOIR -f ${oln}.fasta -v -z 1 -c 10 -r 4 --generate-report --conf ${basedir}/CHOIR.cfg --multiprocess --ignore-templated --allow-monomers
    cd ${basedir}
    #
  fi
done
