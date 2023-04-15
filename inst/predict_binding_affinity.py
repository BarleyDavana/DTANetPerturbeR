import sys
import pandas as pd
from io import StringIO
from Bio import SeqIO
from Bio import Entrez
from time import sleep
from DeepPurpose import DTI as models

def predict_binding_affinity(dti_file, result_folder, fasta_file=None, model='MPNN_CNN_BindingDB_IC50', get_sequence=None, email=None, api_key=None):
    """
    The drug-target interaction file is read, protein sequences are added
    and binding affinity is predicted using a pre-trained DeepPurpose model.

    Args:
        dti_file (str): Path to the data.frame file containing DTIs.
        result_folder (str): The path to the folder where the results are saved.
        fasta_file (str): Path to the FASTA file containing all protein sequences in the DTI file.
        model (str): The DeepPurpose model type used, which defaults to 'MPNN_CNN_BindingDB_IC50'.
        get_sequence (bool): Whether to retrieve protein sequences online (True) or use provided FASTA file (False).
        email (str): For online access to protein sequences, an NCBI email address is required.
        api_key (str): To obtain the protein sequence online, you need to provide your own NCBI API key.

    Returns:
        None.
    """
    # Check if fasta_file is provided
    if fasta_file:
        get_sequence = False
    else:
        get_sequence = True
    
    # Read drug-target interaction files and add protein sequences
    add_sequence_to_dti(dti_file, fasta_file, dti_file, get_sequence)

    # Predicting binding affinity using the DeepPurpose model
    df = pd.read_table(dti_file, sep="\t")
    seq = df['SEQUENCE'].values.tolist()
    smiles = df['SMILES'].values.tolist()
    tar_name = df['GENENAME'].values.tolist()
    drug_name = df['DRUGNAME'].values.tolist()
    net = models.model_pretrained(model=model)
    y_pred = models.virtual_screening(smiles, seq, net, drug_name, tar_name, 
                              result_folder=result_folder, verbose=False)


def add_sequence_to_dti(dti_file, fasta_file, out_file, get_sequence=False, email=None, api_key=None):
    """
    This function reads a data.frame file containing drug-target interactions (DTIs), extracts protein sequences from
    a FASTA file based on the target UniProt IDs, and adds the sequences to the DTI data.frame. The updated 
    data.frame is then saved to a file.

    Parameters:
        dti_file (str): Path to the data.frame file containing DTIs. It must include a "UNIPROTID" column.
        fasta_file (str): Path to the FASTA file containing all protein sequences in the DTI file.
            The file must include sequences for all UniProt IDs in the dti_file.
        out_file (str): Path to the updated data.frame file.
        get_sequence (bool): Whether to obtain protein sequences online. Default is False.
        email (str): Email address to use when accessing the NCBI Entrez database online.
        api_key (str): API key to use when accessing the NCBI Entrez database online.

    Returns:
        None.
    """
    # Read the data.frame file
    df = pd.read_table(dti_file, sep="\t")

    # Read all protein sequences from the FASTA file
    sequences = {}
    with open(fasta_file) as handle:
        for record in SeqIO.parse(handle, "fasta"):
            uniprot_id = record.id.split("|")[1]
            sequences[uniprot_id] = str(record.seq)

    # If online sequence retrieval is requested, call the get_sequence function
    if get_sequence:
        ids = df["UNIPROTID"].tolist()
        sequences = get_sequence(ids, email=email, api_key=api_key)
        sequences_dict = dict(zip(ids, sequences))
    else:
        sequences_dict = sequences

    # Add protein sequences to the data.frame
    sequences_list = []
    for uniprot_id in df["UNIPROTID"]:
        sequence = sequences_dict.get(uniprot_id, "NA")
        sequences_list.append(sequence)

    df["SEQUENCE"] = sequences_list
    
    # Remove rows where the sequence is "NA"
    df = df[df["SEQUENCE"] != "NA"]
    
    # Save the updated data.frame
    df.to_csv(out_file, sep="\t", index=False)


def get_sequence(ids, max_retry=3, retry_delay=2):
    """
    This function retrieves protein sequences from NCBI's Entrez database given a list of UniProt IDs.
    
    Parameters:
        ids (list): A list of UniProt IDs.
        max_retry (int): The maximum number of times to retry a failed request.
        retry_delay (int): The number of seconds to wait between retries.
    
    Returns:
        list: A list of protein sequences corresponding to the input UniProt IDs.
    """
    sequences = []
    # Split the IDs into batches of 100 for batch processing
    for i in range(0, len(ids), 100):
        batch_ids = ids[i:i+100]
        batch_sequences = []
        # Retrieve sequences for each ID in the batch
        for uniprot_id in batch_ids:
            retry_count = 0
            # Retry failed requests up to max_retry times
            while retry_count < max_retry:
                try:
                    # Fetch sequence from Entrez database
                    handle = Entrez.efetch(db="protein", id=uniprot_id, rettype="fasta", retmode="text")
                    record = SeqIO.read(handle, "fasta")
                    handle.close()
                    # Add sequence to batch_sequences list
                    batch_sequences.append(str(record.seq))
                    break
                except Exception as e:
                    retry_count += 1
                    print(f"Error fetching sequence for {uniprot_id}: {e}. Retrying in {retry_delay} seconds...")
                    sleep(retry_delay)
            # If all retries failed, add "NA" to batch_sequences list
            if retry_count >= max_retry:
                batch_sequences.append("NA")
        # Add batch_sequences to sequences list
        sequences += batch_sequences
    # Return list of sequences
    return sequences

if __name__ == '__main__':
    dti_file = sys.argv[1]
    output_folder = sys.argv[2]
    fasta_file = sys.argv[3]
    predict_binding_affinity(dti_file, output_folder, fasta_file)
