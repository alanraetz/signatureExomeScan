# import io
import mhcflurry.fasta

# Define alleles for each sample
alleles={
    "my-sample": [ 'HLA-A*01:01',
'HLA-A*02:01',
'HLA-A*02:03',
'HLA-A*02:06',
'HLA-A*03:01',
'HLA-A*11:01',
'HLA-A*23:01',
'HLA-A*24:02',
'HLA-A*26:01',
'HLA-A*30:01',
'HLA-A*30:02',
'HLA-A*31:01',
'HLA-A*32:01',
'HLA-A*33:01',
'HLA-A*68:01',
'HLA-A*68:02',
'HLA-B*07:02',
'HLA-B*08:01',
'HLA-B*08:01',
'HLA-B*15:01',
'HLA-B*35:01',
'HLA-B*40:01',
'HLA-B*44:02',
'HLA-B*44:03',
'HLA-B*51:01',
'HLA-B*53:01',
'HLA-B*57:01',
'HLA-B*58:01' ],
 # ["A0201", "A0301", "B0702", "C0802"],
}

for suffix in range(97,417):

    # with io.open("file_" + str(i) + ".dat", 'w', encoding='utf-8') as f:

    srcFile = "/content/drive/MyDrive/vaccinePeptideExomeScan/split/exomeScanPeptides_" + str(suffix) + ".fasta"
    
    peptides = mhcflurry.fasta.read_fasta_to_dataframe(srcFile).set_index("sequence_id")

    # resultFile = "/content/drive/MyDrive/vaccinePeptideExomeScan/results/colab_10k_" + str(suffix) + ".csv"
    resultFile = "colab_10k_" + str(suffix) + ".csv"

    print("processing " + resultFile)

    driveResultFile = "/content/drive/MyDrive/vaccinePeptideExomeScan/results/" + resultFile

    # Predict across protein sequences and return peptides with predicted affinity
    # less than 100 nM.
    results2 = predictor.predict_sequences(
        sequences=peptides.sequence.to_dict(),
        alleles=alleles,
        result="filtered",
        comparison_quantity="affinity",
        filter_value=100)
    
    results2.to_csv(driveResultFile)
    # files.download(results2.to_csv())


'HLA-A*01:01',
'HLA-A*02:01',
'HLA-A*02:03',
'HLA-A*02:06',
'HLA-A*03:01',
'HLA-A*11:01',
'HLA-A*23:01',
'HLA-A*24:02',
'HLA-A*26:01',
'HLA-A*30:01',
'HLA-A*30:02',
'HLA-A*31:01',
'HLA-A*32:01',
'HLA-A*33:01',
'HLA-A*68:01',
'HLA-A*68:02',
'HLA-B*07:02',
'HLA-B*08:01',
'HLA-B*08:01',
'HLA-B*15:01',
'HLA-B*35:01',
'HLA-B*40:01',
'HLA-B*44:02',
'HLA-B*44:03',
'HLA-B*51:01',
'HLA-B*53:01',
'HLA-B*57:01',
'HLA-B*58:01'