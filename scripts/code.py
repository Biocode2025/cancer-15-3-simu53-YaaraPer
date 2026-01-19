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
def Comp_seq(old,new):
  '''
  הפונקציה בודקת כמה הבדלים קיימים בין הרצפים השונים ומחזירה את מספר ההבדלים.
  מקבלת: old,new.
  מחזירה: num_differences.
  '''
  num_differences = 0
  
  for i in range(min(len(old), len(new))):
    if old[i] != new[i]:
      num_differences = num_differences + 1
      
  return num_differences  
#------------------------------------------------

def Insert_DNA(seq):
  '''
  הפונקציה תכניס במיקום אקראי לרצף ה- DNA של הגן נוקלאוטיד נוסף.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  nucleotide_list = ['T','G','C','A']
  
  rand_nucleotide = random.choice(nucleotide_list)
  rand_num = random.randrange(0,len(seq))
  
  change_genome = seq[0:rand_num]+ rand_nucleotide + seq[rand_num:]
  
  return change_genome
#------------------------------------------------
  
def Delete_DNA(seq):
  '''
  הפונקציה תחסיר נוקלאוטיד במיקום רנדומאלי.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  rand_num = random.randrange(0,len(seq))
  rand_nucleotide = seq[rand_num]
  
  change_genome = seq[:rand_num] + seq[rand_num + 1:]

  return change_genome
#------------------------------------------------

def Mutate_DNA(seq):
  '''
  הפונקציה בוחרת מיקום אקראי ברצף גנום ה- HIV ומחליפה במיקום זה את הנוקלאוטיד באופן רנדומלי ל- A/C/G/T.
  מקבלת: seq.
  מחזירה:change_genome
  '''
  nucleotide_list = ['T','G','C','A']
  
  rand_nucleotide = random.choice(nucleotide_list)
  rand_num = random.randrange(0,len(seq))
  
  if seq[rand_num] != rand_nucleotide:
    change_genome = seq[0:rand_num]+ rand_nucleotide + seq[(rand_num+1):]
  
  else:
    nucleotide_list.remove(rand_nucleotide)
    rand_nucleotide = random.choice(nucleotide_list)
    change_genome = seq[0:rand_num]+ rand_nucleotide + seq[(rand_num+1):]
  return change_genome 

# תוכנית ראשית.
global RNA_codon_table
RNA_codon_table = {}

num_gen = 1000

# פתיחת הקבצים
p53_seq = open('data/human_p53_coding.txt', 'r')
codon_file = open('data/codon_AA (1).txt', 'r')

# קריאה לפונקציה
Read_dict(codon_file)

p53_genome = ""
new_HIV_genome = ""
# קריאת הקובץ
for line in p53_seq:
  line = line.rstrip('\r\n')
  if line == "":
    continue
  # רצף ה DNA מופיע בשורות שאינן מתחילות בסימן "<" לכן "נדלג" על שורה זו
  if line[0] == ">":
    continue
  
  p53_genome = p53_genome + line

# קריאה לפונקציות- שעתוק ותרגום הרצף.
gene_as_RNA = DNA_RNA_Cod(p53_genome)
old_protein = RNA_prot(gene_as_RNA)

is_changed = True

while (is_changed)
  num = random.randrange(1,101)

  # מוטציה של החלפת בסיס
  if num <= 98:
    p53_genome = Mutate_DNA(p53_genome)   # האם להחשיב כמוטציה גם אם היא שקטה?

  # מוטציה של הוספת בסיס עד שלושה בסיסים
  elif num == 1:
    num_bases = random.randrange(1,4)
    if num_bases == 1:
      p53_genome = Insert_DNA(p53_genome)
    elif num_bases == 2:
      for i in range(num_bases):
        p53_genome = Insert_DNA(p53_genome)
    else:
      for i in range(num_bases):
        p53_genome = Insert_DNA(p53_genome)
  
  # מוטציה של הוספת בסיס עד שלושה בסיסים
  else:
    num_bases = random.randrange(1,4)
    if num_bases == 1:
      p53_genome = Delete_DNA(p53_genome)
    elif num_bases == 2:
      for i in range(num_bases):
        p53_genome = Delete_DNA(p53_genome)
    else:
      for i in range(num_bases):
        p53_genome = Delete_DNA(p53_genome)
  
  # קריאה לפונקציות- שעתוק ותרגום הרצף.
  gene_as_RNA = DNA_RNA_Cod(p53_genome)
  new_protein = RNA_prot(gene_as_RNA)

  if Comp_seq(old_protein, new_protein) > 0:
    is_changed = False
  
  else:
    num_iteration = num_iteration + 1


    






# התוכנית מבצעת שלוש מוטציות נקודתיות בצורה אקראית לאורך הרצף.
for i in range(3):
  num = random.randrange(3)
  if num == 1:
    p53_genome = Mutate_DNA(p53_genome)
  elif num == 2:
    p53_genome = Insert_DNA(p53_genome)
  else:
    p53_genome = Delete_DNA(p53_genome)

# קריאה לפונקציות- שעתוק ותרגום הרצף.
gene_as_RNA = DNA_RNA_Cod(p53_genome)
new_protein = RNA_prot(gene_as_RNA)

