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

def Insert_DNA(seq, num):
  '''
  הפונקציה תכניס במיקום אקראי לרצף ה- DNA של הגן נוקלאוטיד/ים נוסף.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  nucleotide_list = ['T','G','C','A']
  change_genome = ""


  if num != 0 and num != 1:
    for i in range(num):
      change_genome = change_genome + random.choice(nucleotide_list)


  if len(seq) == 0:
     return seq


  rand_nucleotide = random.choice(nucleotide_list)
  rand_num = random.randrange(0,len(seq))
 
  change_genome = seq[0:rand_num]+ rand_nucleotide + seq[rand_num:]
 
  return change_genome
#------------------------------------------------
 
def Delete_DNA(seq, num):
  '''
  הפונקציה תחסיר נוקלאוטיד/ים במיקום רנדומאלי.
  מקבלת: seq.
  מחזירה: change_genome.
  '''
  if len(seq) == 0 and num < len(seq):
    return seq


  rand_num = random.randrange(0,len(seq))
  rand_nucleotide = seq[rand_num]
 
  change_genome = seq[:rand_num] + seq[rand_num + num:]
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
#------------------------------------------------

# תוכנית ראשית.
global RNA_codon_table
RNA_codon_table = {}

# הגדרת משתנים
num_gen = 1000
avg = 0
p53_genome = ""
num_iteration = 0
num_mutate = 0
# הגדרת רשימה
iteration_list = []

# פתיחת הקבצים
p53_seq = open('data/human_p53_coding (1).txt', 'r')
codon_file = open('data/codon_AA (1).txt', 'r')
# קריאה לפונקציה
Read_dict(codon_file)

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
old_protein = RNA_prot(DNA_RNA_Cod(p53_genome))

# לולאה חיצונית- עובדת לפי מספר הדורות
for h in range(num_gen):
  is_changed = True

  # לולאה פנימית- מדמה את תהליך התרגום של החלבון, בה מתרחשות המוטציות והבדיקה של כמות המוטציות בכל דור.
  while (is_changed):
    num_iteration = num_iteration + 1
    
    # הגרלת מספרים
    Mutate_rnd_num = random.randint(1,100)
    frequency_rnd_num = random.randint(1,10000)
    
    # התדירות האפקטיבית להתרחשות מוטציה היא  1 ל 10000 אירועי תרגום לכן אם נקבל בהגרלה מספר אחד מבין 10000 מספרים תתרחש מוטציה כלשהי
    if frequency_rnd_num == 1:
      # מוטציה של החלפת בסיס
      if Mutate_rnd_num <= 98:
        p53_genome = Mutate_DNA(p53_genome)  

      # מוטציה של הוספת בסיס עד שלושה בסיסים
      elif Mutate_rnd_num == 99:
        num_bases = random.randrange(1,4)
        p53_genome = Insert_DNA(p53_genome, num_bases)
      
      # מוטציה של הוספת בסיס עד שלושה בסיסים
      else:
        num_bases = random.randrange(1,4)
        p53_genome = Delete_DNA(p53_genome, num_bases)
    
    # קריאה לפונקציות- שעתוק ותרגום הרצף.
    new_protein = RNA_prot(DNA_RNA_Cod(p53_genome))

    num_mutate = num_mutate + Comp_seq(old_protein, new_protein)
    # היות ומספיקה מוטציה אחת לעצירת הדור
    if num_mutate >= 1:
      is_changed = False
      num_mutate = 0


  # קריאה לפונקציות- שעתוק ותרגום הרצף.      
  old_protein = new_protein
  # סכימת מספר האיטרציות שלקח ללולאה הפנימית לעשות עד שנוצרה מוטציה (לא שקטה).
  iteration_list.append(num_iteration)
  num_iteration = 0

# חישוב סכום מספר הדורות לקבלת מוטציה בגן (שאינה שקטה).
total = sum(iteration_list)

# חישוב ממוצע הדורות לקבלת מוטציה בגן (שאינה שקטה).
avg = total / num_gen
print(avg)
