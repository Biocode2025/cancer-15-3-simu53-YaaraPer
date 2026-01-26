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
