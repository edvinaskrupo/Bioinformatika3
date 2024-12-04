import matplotlib.pyplot as plt
import numpy as np

def cg_distribution_from_fastq(fastq_file):
    cg_fractions = []
    
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        
        # Pradėkime nuo pirmo read'o ir pereikime per visus
        for i in range(0, len(lines), 4):  # FASTQ įrašas turi keturias eilutes
            sequence = lines[i + 1].strip()  # Sekos eilutė
            cg_count = sequence.count('C') + sequence.count('G')  # C ir G nukleotidai
            cg_fraction = cg_count / len(sequence)  # C/G dalis
            cg_fractions.append(cg_fraction)
    
    return cg_fractions

fastq_file = 'reads_for_analysis.fastq'
cg_fractions = cg_distribution_from_fastq(fastq_file)

# Sukuriame grafiką
plt.hist(cg_fractions, bins=50)
plt.xlabel('C/G nukleotidų dalis (procentais)')
plt.ylabel('Read’ų skaičius')
plt.title('C/G nukleotidų pasiskirstymas read’uose')
plt.show()

