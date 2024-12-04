def detect_quality_encoding(quality_scores):
    # Paimame pirmąją kokybės reikšmę
    score = ord(quality_scores[0])
    
    # Tikriname pagal Phred+33 ar Phred+64 kodavimus
    if score >= 33 and score <= 73:  # Phred+33 (Sanger, Illumina 1.8+)
        return 'Phred+33 (Sanger, Illumina 1.8+)'
    elif score >= 64 and score <= 104:  # Phred+64 (Illumina 1.3+, 1.5+)
        return 'Phred+64 (Illumina 1.3+, 1.5+)'
    elif score >= 59 and score <= 104:  # Solexa+64
        return 'Solexa+64'
    else:
        return 'Nežinomas kodavimas'

def determine_quality_encoding_from_fastq(fastq_file):
    with open(fastq_file, 'r') as f:
        lines = f.readlines()
        
        # Kokybės eilutė yra kas ketvirta eilutė (4-tas simbolis)
        quality_scores = lines[3].strip()  # Paimame ketvirtą eilutę (kokybės reikšmes)
        
        encoding = detect_quality_encoding(quality_scores)
        return encoding

fastq_file = 'reads_for_analysis.fastq'
encoding = determine_quality_encoding_from_fastq(fastq_file)
print(f"Naudojamas kokybės kodavimas: {encoding}")
