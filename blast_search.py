from Bio import SeqIO
from Bio.Blast import NCBIWWW
from Bio.Blast import NCBIXML
import pandas as pd
import datetime

# parametrai
file_path = "reads_for_analysis.fastq"
distribution_peaks = [0.35, 0.55, 0.7] # pikai
nr_of_reads_per_peak = 5  # sekų skaičius per piką
truncate_length = 50 # sutrumpiname sekas greitesnei BLAST paieškai

def log_message(message):
    print(f"[{datetime.datetime.now()}] {message}")

def get_reads_from_peaks(file_path, peaks, reads_per_peak):
    selected_reads = []
    peak_reads_count = {peak: 0 for peak in peaks}  # skaičiuojame kiek sekų nuskaitėme nuo piko

    for i, record in enumerate(SeqIO.parse(file_path, "fastq")):
        if len(selected_reads) >= len(peaks) * reads_per_peak:
            break
        
        for peak in peaks:
            if i % 100 == int(peak * 100) and peak_reads_count[peak] < reads_per_peak:
                selected_reads.append(record)
                peak_reads_count[peak] += 1
                break

    return selected_reads

log_message("Starting sequence selection...")
reads_for_analysis = get_reads_from_peaks(file_path, distribution_peaks, nr_of_reads_per_peak)

# atliekame BLAST paiešką
results = []
for read in reads_for_analysis:
    truncated_sequence = read.seq[:truncate_length]
    log_message(f"Starting BLAST search for read: {read.id} (first {truncate_length} bp)")
    result_handle = NCBIWWW.qblast("blastn", "nr", truncated_sequence, entrez_query="bacteria[organism]")
    blast_record = NCBIXML.read(result_handle)
    
    if blast_record.alignments:
        best_hit = blast_record.alignments[0]
        organism_name = best_hit.hit_def.split("[")[-1].strip("]")
        results.append({"Read ID": read.id, "Organism": organism_name})
        log_message(f"Found organism: {organism_name}")
    else:
        log_message(f"No matches found for read: {read.id}")
        results.append({"Read ID": read.id, "Organism": "No match"})

# išvedame rezultatus
df = pd.DataFrame(results)
df.to_csv("blast_results.csv", index=False)
log_message("BLAST search completed. Results saved to 'blast_results.csv'.")
