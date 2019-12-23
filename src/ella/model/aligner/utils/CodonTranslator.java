package ella.model.aligner.utils;

public class CodonTranslator {

	public static char translateCodon(String codon) {

		if (codon.length() != 3)
			return 'X';

		int a = getAccess(codon.charAt(0));
		int b = getAccess(codon.charAt(1));
		int c = getAccess(codon.charAt(2));

		return codons[a][b][c];

	}

	public static String translate(String dna, int frame) {
		StringBuilder prot = new StringBuilder(dna.length() / 3);
		StringBuilder codon = new StringBuilder(3);
		for (int i = Math.abs(frame) - 1; i < dna.length(); i++) {
			codon.append(dna.charAt(i));
			if (codon.length() == 3) {
				char p = translateCodon(codon.toString());
				prot = prot.append(p);
				codon = new StringBuilder(3);
			}
		}
		return prot.toString();
	}

	public static char translateCodon(char c1, char c2, char c3) {
		int a = getAccess(c1);
		int b = getAccess(c2);
		int c = getAccess(c3);
		return codons[a][b][c];
	}

	private static char[][][] codons;

	static {

		codons = new char[26][26][26];
		for (int a = 0; a < 26; a++) {
			for (int b = 0; b < 26; b++) {
				for (int c = 0; c < 26; c++) {
					codons[a][b][c] = 'X';
				}
			}
		}

		addCodon("GGT", 'G');
		addCodon("GGU", 'G');
		addCodon("GGC", 'G');
		addCodon("GGA", 'G');
		addCodon("GGG", 'G');

		addCodon("GCT", 'A');
		addCodon("GCU", 'A');
		addCodon("GCC", 'A');
		addCodon("GCA", 'A');
		addCodon("GCG", 'A');

		addCodon("GTT", 'V');
		addCodon("GUU", 'V');
		addCodon("GTC", 'V');
		addCodon("GUC", 'V');
		addCodon("GTA", 'V');
		addCodon("GUA", 'V');
		addCodon("GTG", 'V');
		addCodon("GUG", 'V');

		addCodon("TTA", 'L');
		addCodon("UUA", 'L');
		addCodon("TTG", 'L');
		addCodon("UUG", 'L');
		addCodon("CTT", 'L');
		addCodon("CUU", 'L');
		addCodon("CTC", 'L');
		addCodon("CUC", 'L');
		addCodon("CTA", 'L');
		addCodon("CUA", 'L');
		addCodon("CTG", 'L');
		addCodon("CUG", 'L');

		addCodon("ATT", 'I');
		addCodon("AUU", 'I');
		addCodon("ATC", 'I');
		addCodon("AUC", 'I');
		addCodon("ATA", 'I');
		addCodon("AUA", 'I');

		addCodon("TCT", 'S');
		addCodon("UCU", 'S');
		addCodon("TCC", 'S');
		addCodon("UCC", 'S');
		addCodon("TCA", 'S');
		addCodon("UCA", 'S');
		addCodon("TCG", 'S');
		addCodon("UCG", 'S');
		addCodon("AGT", 'S');
		addCodon("AGU", 'S');
		addCodon("AGC", 'S');

		addCodon("ACT", 'T');
		addCodon("ACU", 'T');
		addCodon("ACC", 'T');
		addCodon("ACA", 'T');
		addCodon("ACG", 'T');

		addCodon("GAT", 'D');
		addCodon("GAU", 'D');
		addCodon("GAC", 'D');

		addCodon("GAA", 'E');
		addCodon("GAG", 'E');

		addCodon("AAT", 'N');
		addCodon("AAU", 'N');
		addCodon("AAC", 'N');

		addCodon("CAA", 'Q');
		addCodon("CAG", 'Q');

		addCodon("AAG", 'K');
		addCodon("AAA", 'K');

		addCodon("CGT", 'R');
		addCodon("CGU", 'R');
		addCodon("CGC", 'R');
		addCodon("CGA", 'R');
		addCodon("CGG", 'R');
		addCodon("AGA", 'R');
		addCodon("AGG", 'R');

		addCodon("CAT", 'H');
		addCodon("CAU", 'H');
		addCodon("CAC", 'H');

		addCodon("TTT", 'F');
		addCodon("UUU", 'F');
		addCodon("TTC", 'F');
		addCodon("UUC", 'F');

		addCodon("TGT", 'C');
		addCodon("UGU", 'C');
		addCodon("TGC", 'C');
		addCodon("UGC", 'C');

		addCodon("TGG", 'W');
		addCodon("UGG", 'W');

		addCodon("TAT", 'Y');
		addCodon("UAU", 'Y');
		addCodon("TAC", 'Y');
		addCodon("UAC", 'Y');

		addCodon("ATG", 'M');
		addCodon("AUG", 'M');

		addCodon("CCT", 'P');
		addCodon("CCU", 'P');
		addCodon("CCC", 'P');
		addCodon("CCA", 'P');
		addCodon("CCG", 'P');

	}

	private static void addCodon(String codon, char r) {
		int a = getAccess(codon.charAt(0));
		int b = getAccess(codon.charAt(1));
		int c = getAccess(codon.charAt(2));
		codons[a][b][c] = r;
	}

	private static int getAccess(char c) {
		return c > 122 ? c - 32 - 65 : c - 65;
	}

}
