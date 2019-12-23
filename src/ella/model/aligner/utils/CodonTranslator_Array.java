package ella.model.aligner.utils;

public class CodonTranslator_Array {

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

		addCodon("GGU", 'G');
		addCodon("GGT", 'G');
		addCodon("GGC", 'G');
		addCodon("GGA", 'G');
		addCodon("GGG", 'G');

		addCodon("GCU", 'A');
		addCodon("GCT", 'A');
		addCodon("GCC", 'A');
		addCodon("GCA", 'A');
		addCodon("GCG", 'A');

		addCodon("GUU", 'V');
		addCodon("GTT", 'V');
		addCodon("GUC", 'V');
		addCodon("GTC", 'V');
		addCodon("GUA", 'V');
		addCodon("GTA", 'V');
		addCodon("GUG", 'V');
		addCodon("GTG", 'V');

		addCodon("UUA", 'L');
		addCodon("TTA", 'L');
		addCodon("UUG", 'L');
		addCodon("TTG", 'L');
		addCodon("CUU", 'L');
		addCodon("CTT", 'L');
		addCodon("CUC", 'L');
		addCodon("CTC", 'L');
		addCodon("CUA", 'L');
		addCodon("CTA", 'L');
		addCodon("CUG", 'L');
		addCodon("CTG", 'L');

		addCodon("AUU", 'I');
		addCodon("ATT", 'I');
		addCodon("AUC", 'I');
		addCodon("ATC", 'I');
		addCodon("AUA", 'I');
		addCodon("ATA", 'I');

		addCodon("UCU", 'S');
		addCodon("TCT", 'S');
		addCodon("UCC", 'S');
		addCodon("TCC", 'S');
		addCodon("UCA", 'S');
		addCodon("TCA", 'S');
		addCodon("UCG", 'S');
		addCodon("TCG", 'S');
		addCodon("AGU", 'S');
		addCodon("AGT", 'S');
		addCodon("AGC", 'S');

		addCodon("ACU", 'T');
		addCodon("ACT", 'T');
		addCodon("ACC", 'T');
		addCodon("ACA", 'T');
		addCodon("ACG", 'T');

		addCodon("GAU", 'D');
		addCodon("GAT", 'D');
		addCodon("GAC", 'D');

		addCodon("GAA", 'E');
		addCodon("GAG", 'E');

		addCodon("AAU", 'N');
		addCodon("AAT", 'N');
		addCodon("AAC", 'N');

		addCodon("CAA", 'Q');
		addCodon("CAG", 'Q');

		addCodon("AAG", 'K');
		addCodon("AAA", 'K');

		addCodon("CGU", 'R');
		addCodon("CGT", 'R');
		addCodon("CGC", 'R');
		addCodon("CGA", 'R');
		addCodon("CGG", 'R');
		addCodon("AGA", 'R');
		addCodon("AGG", 'R');

		addCodon("CAU", 'H');
		addCodon("CAT", 'H');
		addCodon("CAC", 'H');

		addCodon("UUU", 'F');
		addCodon("TTT", 'F');
		addCodon("UUC", 'F');
		addCodon("TTC", 'F');

		addCodon("UGU", 'C');
		addCodon("TGT", 'C');
		addCodon("UGC", 'C');
		addCodon("TGC", 'C');

		addCodon("UGG", 'W');
		addCodon("TGG", 'W');

		addCodon("UAU", 'Y');
		addCodon("TAT", 'Y');
		addCodon("UAC", 'Y');
		addCodon("TAC", 'Y');

		addCodon("AUG", 'M');
		addCodon("ATG", 'M');

		addCodon("CCU", 'P');
		addCodon("CCT", 'P');
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
		return c >= 97 ? c - 32 - 65 : c - 65;
	}

	public static char translateCodon(String codon) {
		if (codon.length() != 3)
			return 'X';
		int a = getAccess(codon.charAt(0));
		int b = getAccess(codon.charAt(1));
		int c = getAccess(codon.charAt(2));
		return codons[a][b][c];

	}

	public static String translate(String dna, int frame) {
//		String rna = dna.replaceAll("T", "U");
		StringBuffer prot = new StringBuffer("");
		String codon = "";
		for (int i = Math.abs(frame) - 1; i < dna.length(); i++) {
			codon = codon.concat("" + dna.charAt(i));
			if (codon.length() == 3) {
				char p = translateCodon(codon);
				prot = prot.append(p);
				codon = "";
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

}
