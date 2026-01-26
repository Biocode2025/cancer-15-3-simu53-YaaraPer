import random

def DNA_RNA_Cod(seq):
  '''
  הפונקציה דואגת שהאותיות תהיינה אחידות (אותיות גדולות) והופכת את רצף ה- DNA המקודד לרצף RNA.
  מקבלת: seq.
  מחזירה: RNA_seq.
  '''
  RNA_seq=""


  for line in seq:
    line = line.upper()
    line = line.rstrip('\r\n')
   
    if line[0] != ">":
      rna_line = line.replace("T","U")
      RNA_seq += rna_line
   
  return RNA_seq
#------------------------------------------------
 
def Read_dict(fl):
  '''
  הפונקציה קוראת לתוךdictionary  את המיפוי בין הקודונים לחומצות אמינו מהקובץ.
  מקבלת: fl.
  מחזירה: RNA_codon_table.
  '''
  global RNA_codon_table


  for line in fl:
    line = line.rstrip('\r\n')
    line_list = line.split()
   
    codon_key = line_list[0]
    amino_acids_value = line_list[1]
   
    RNA_codon_table[codon_key] = amino_acids_value
   
  return RNA_codon_table
#------------------------------------------------
 
def RNA_prot(seq):
  '''
  הפונקציה מתרגמת את רצף ה- RNA לרצף חלבון.
  מקבלת: seq.
  מחזירה: protein_seq.
  '''
  protein_seq = ""
  start = False
  for i in range(0, len(seq), 3):
    codon = seq[i:i+3]
    if len(codon) < 3:
      break
   
    curr_amino_acid = RNA_codon_table.get(codon,"-")
   
    if curr_amino_acid == "M":
      start = True
     
    if start:
      if curr_amino_acid == "*":
        protein_seq = protein_seq + "*"
        break
     
      protein_seq += curr_amino_acid
 
  return protein_seq
#------------------------------------------------
