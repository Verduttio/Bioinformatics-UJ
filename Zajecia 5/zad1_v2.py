from Bio import pairwise2
from Bio import SeqIO
from Bio.Align import substitution_matrices


if __name__ == '__main__':
    subs_mat = substitution_matrices.load('BLOSUM62')

    seq11 = "PLEASANTLY"
    seq12 = "MEANLY"

    seq1 = "SWFQIMVDLSLFWFHPYLHEKNVDAKTQDATMSQTCYAHSRGETRRRRCTKHPDCPHLPTMAFQEPVQPAGLHWGGVHLAKEPLCNSQQHQKIFQTSIKQPFVSGQVHDNGNTMVIECGEESHAWPPHRTIHSFHIDVCQWIKEMTSQFVYKNSEHLWEESPENGENLFILHCCSAEHKMQRPSSNHNFFAITREIVEDRTLGPLQASAGIWGIVIPLLPRCQLMMNWNHYDWWAWWAMIPAYGGGETLRPHSMRNCLEDGIKCNHVCFSNLIKFGHWDWRYNKKKHPGEWNSRYMLHRAVAETNIVKRHMFCFHLIETKCVALHQLLELTKNFLFAMFVDDVTPLWDACHECCIKWKAIYMMLYANQSEKLVYKKMSVQVVKNWFPPDVCSGCKSWRACEYMLYFMVNELTYMTCRMHVNPVCVRHWFLVKCSNFCAWFLQRYNHRKAHKGPHGVTQLNMKYRALKAEWVNANWYKMKMTSMVKLVTQSHFMYQTDFFIIMHSIKGAWGYWGYNNTHNIVGEDYNPSRNQWLPMDPQKGWLQNFMCCCWMYGSMIINFPKILPCYPCPPRDHVLKFFINEMFGFMDMIFKCSNVKNDHLLACQLITPPGPKVGLYHLIVYVMNSYINVQKWEFCFTKREEGTLFCLDWEEQHYEYIRMNNKVMYSIWWAVNTNYGPIVCVTWAIYTWRMRADHNENFFGNNGNNMMHTEQKVEHGEAFVIHRCQNCCWGRGPKEIMVAMLPVGMMPGGQHEQSCQTEKCGITIRTNHHVFCKGVKAMEAEMGNCYARVNGYMGHSEYNSMNETQKPDAWFVMQGFESPMIYESFYTY"
    seq2 = "SWDQIMVHPYKTMDATSSQTCCEVATIDGAHSRGETRRRHPDCPSWFQMPAGLHWGGVHL AKEPLCNSQQHQKIIKQPFVSGQVHDNGNTMVFEAGEESHPPHSCPFHGDVEMKSQFVKDHDWMQIAKNSEHLWPENYCVNMGASAAYPRTYPHDSSSSHNFFAITREIVEVRWLGPLQAMVQGIVIPLLPRCKMGSPLMMNYFPSSWQLLWAMINAYDGGGGIPQRETLFPHSMRNCLEVSSHYGIKCNHVYFIKFIHWDWRYNKKKSPWNSRPQMWFYKMGDNNNLLHVAVAETSIVKRHMFCFHLIECLDKFYKCWFTPWKHFQKNFLFWDACHECPKKWKVIYMMLYANQSEKFVYKKKSVQVVKNWQYSECPRVDNTNDNGPCLCSGCKSWRACEYMLYFCTYMTCAMNCVSPNSPGCTGWGDPCIFKWFLQRYNHRKAHKGPHRVTQLNMKYLKAEWVGLRWNWYKGNKMTSMVMYQPRMFDFFIINFAWCLTHSIKGAWGFWMYNEHFKIEEKTHAKWKDHSTKCKNEDDTNPSRNVWYPMCPQKGVGQNFMCCWMYGSMCINFPKILPCPPRDHVLKTFINEMFGFMDAIFKCSNYKNDHLLGHQDSWLYHVIVYVMNSYINVQKWEFAFTKREEGTLHCWPHQGSHCDWEESHYEYCPMNNKKMYSIWWGVNGNYGPVICWTWAIYTCRMRADQQWCPGNENFFGNNGNDMMHTEQKVEHGCKTYSMAFVIHRCQNCNWGRGPKEIMVAMLPPMRLECQTEKCGITWRTNHHKEAEMGNCYARVTVYMGHSEYNSDNETQKPCAWHVMQGFESFNMIYESFYTY"
    alignments = pairwise2.align.globalds(seq1, seq2, subs_mat, -5, -5)
    len(alignments)
    print(pairwise2.format_alignment(*alignments[0]))