import os
import random
import time
import jpredapi
import shutil


def generate_random_sequences(seq, n):
    random_sequences = []
    for i in range(n):
        permutation = generate_permutation(seq)
        random_sequences.append(permutation)
    return random_sequences


def generate_permutation(seq):
    seq_list = list(seq)
    random.shuffle(seq_list)
    return ''.join(seq_list)


def get_job_id(jp_submit_response):
    return (str(jp_submit_response.content).split("?"))[1].split("\"")[0]


def save_results_to_file(job_id):
    # remove directory if exists
    shutil.rmtree('jpred_sspred/results', ignore_errors=True)
    while not os.path.exists("jpred_sspred/results"):
        jpredapi.get_results(jobid=job_id, results_dir_path="jpred_sspred/results", extract=True)
        time.sleep(10)


def predicted_secondary_structure(job_id):
    with open(f"jpred_sspred/results/{job_id}/{job_id}.jnet", "r") as f:
        file = f.readlines()
        seq = file[1].split(':')[1]
        return seq.replace(',', '')


def make_secondary_structure_prediction(seq):
    result = jpredapi.submit(mode="single", user_format="raw", seq=seq)
    job_id = get_job_id(result)
    save_results_to_file(job_id)
    return predicted_secondary_structure(job_id)


def do_hypothesis(seq, sec_seq):
    hypothesis_results = []
    real_seq_el_count = count_secondary_structure_elements(sec_seq)
    random_sequences = generate_random_sequences(seq, 5)
    for random_sequence in random_sequences:
        predicted_sec_str = make_secondary_structure_prediction(random_sequence)
        print(f"Real sequence: {random_sequence}")
        print(f"Predicted secondary structure: {predicted_sec_str}")
        sec_str_elements = count_secondary_structure_elements(predicted_sec_str)
        print(f"Predicted secondary structure elements count: {sec_str_elements}")
        hypothesis_result = hypothesis_term(real_seq_el_count, sec_str_elements)
        print(f"Real sequence secondary structure elements count: {real_seq_el_count}")
        print(f"Hypothesis result (predicted_count < real_count*0.8: {hypothesis_result}")
        hypothesis_results.append(hypothesis_result)
    return hypothesis_results


def count_secondary_structure_elements(seq):
    return len(seq) - seq.count('-')


def hypothesis_term(real_seq_el_count, predicted_seq_el_count):
    return predicted_seq_el_count <= real_seq_el_count * 0.8


def check_hypothesis(hypothesis_results):
    positive = 0
    for result in hypothesis_results:
        if result:
            positive += 1
    negative = len(hypothesis_results) - positive
    print("----------------------------------")
    print("Hypothesis is true if predicted secondary structure elements count is less than or equal 80% of real "
          "sequence secondary structure elements count")
    print(f"Positive results: {positive}")
    print(f"Negative results: {negative}")
    print(f"Positive results percentage: {positive / len(hypothesis_results)}")
    print(f"Is positive percentage higher than 80%? {positive / len(hypothesis_results) > 0.8}")
    print(f"Hypothesis result: {positive / len(hypothesis_results) > 0.8}")


if __name__ == '__main__':
    ubiquitin = "MQIFVKTLTGKTITLEVEPSDTIENVKAKIQDKEGIPPDQQRLIFAGKQLEDGRTLSDYNIQKESTLHLVLRLRGG"
    ubiquitin_secondary_structure = "EEEEEEE--EEEEEEEE-----HHHHHHHHHHHH-----EEEEEE--EEE-----HHHH----EEEEEEEEE----"

    seq_3LIY = "PVIPLDPARRPVIKAQVDTQTSHPKTIEALLDTGADMTVIPIALFSSNTPLKNTSVLGAGGQTQDHFKLTSLPVLIRLPFRTTPIVLTSCLVDTKNNWAIIGRDALQQCQGVLYLP"
    sec_str_3LIY = "-EEE-------EEEEEEE------EEEEEEE-------EEE-HHH-----EEE--EEE--EEE---EEEE---EEEE-------EEE---EEE------EEHHHHHHHH--EEE--"

    hypothesis_results = []
    hypothesis_results.append(do_hypothesis(ubiquitin, ubiquitin_secondary_structure))
    hypothesis_results.append(do_hypothesis(seq_3LIY, sec_str_3LIY))

    check_hypothesis(hypothesis_results)
