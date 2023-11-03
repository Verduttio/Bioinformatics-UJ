def dna_strings_to_array(input_file):
    dnas = []
    dna = ""
    input_file.readline()
    while True:
        line = input_file.readline()
        if not line: #EOF
            dnas.append(dna)
            break
        if not new_dna(line):
            if line[len(line)-1] == "\n":
                line = line[:len(line)-1]
            dna += line

        else:
            dnas.append(dna)
            dna = ""
    return dnas


def new_dna(line):
    return line[0] == ">"


def common_substring(dnas):
    # https://www.geeksforgeeks.org/longest-common-substring-array-strings/
    n = len(dnas)
    s = dnas[0]
    ll = len(s)

    res = ""

    for i in range(ll):
        for j in range(i + 1, ll + 1):
            stem = s[i:j]
            k = 1
            for k in range(1, n):
                if stem not in dnas[k]:
                    break

            if k + 1 == n and len(res) < len(stem):
                res = stem

    return res


if __name__ == '__main__':
    file = open("input.txt", "r")
    dnas = dna_strings_to_array(file)
    print(dnas)
    print(common_substring(dnas))
