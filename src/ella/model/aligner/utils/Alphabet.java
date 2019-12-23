package ella.model.aligner.utils;

public class Alphabet {

	public final static String AA_STRING = "ARNDCQEGHILKMFPSTWYVBJZX*/\\";
	private final static int[] reductionGroupSizes = new int[256];
	private final static char[] reductionArray = new char[256];

	private static String[] AA_GROUPS = { "K", "R", "E", "D", "Q", "N", "C", "G", "H", "I", "L", "V", "M", "F", "Y", "W", "P", "S", "T", "A" };
	// private static String[] AA_GROUPS = { "KREDQN", "C", "G", "H", "ILV", "M", "F", "Y", "W", "P", "STA" };
	// private static String[] AA_GROUPS = { "ILMV", "WYF", "P", "C", "GA", "STQN", "DE", "KRH" };
	// private final static String[] AA_GROUPS = { "ILMV", "FWY", "A", "C", "G", "H", "P", "KR", "ST", "DENQ" };

	static {
		for (int i = 0; i < 256; i++) {
			reductionArray[i] = '*';
			reductionGroupSizes[i] = 1;
		}
		for (int i = 0; i < AA_STRING.length(); i++)
			reductionArray[AA_STRING.charAt(i)] = AA_STRING.charAt(i);
		for (String s : AA_GROUPS) {
			for (int i = 0; i < s.length(); i++) {
				reductionArray[s.charAt(i)] = s.charAt(0);
				reductionGroupSizes[s.charAt(i)] = s.length();
			}
		}
	}

	private static int[] indexArray = new int[256];
	static {
		for (int i = 0; i < 256; i++)
			indexArray[i] = AA_STRING.indexOf('*');
		for (int i = 0; i < AA_STRING.length(); i++)
			indexArray[AA_STRING.charAt(i)] = i;
	}

	public static String getAaString() {
		return AA_STRING;
	}

	public static char reduceCharacter(char c) {
		return reductionArray[c];
	}

	public static int reducePosition(int i) {
		return getIndex(getReducedCharacter(i));
	}

	public static char getReducedCharacter(int i) {
		return reductionArray[getCharacter(i)];
	}

	public static char getCharacter(int i) {
		if (i < AA_STRING.length())
			return AA_STRING.charAt(i);
		return '*';
	}

	public static int getIndex(char c) {
		return indexArray[c];
	}

	public static int getReductionGroupSize(char c) {
		return reductionGroupSizes[c];
	}

	public static String reduceSequence(String s) {
		StringBuilder builder = new StringBuilder(s.length());
		for (int i = 0; i < s.length(); i++)
			builder.append(reduceCharacter(s.charAt(i)));
		return builder.toString();
	}

	public static String getReducedAminoacids() {
		String reduced_aminoacids = "BJZX*";
		for (String group : AA_GROUPS)
			reduced_aminoacids = reduced_aminoacids.concat(group.charAt(0) + "");
		String s = "";
		for (int i = 0; i < AA_STRING.length(); i++) {
			char c = AA_STRING.charAt(i);
			if (reduced_aminoacids.contains(c + ""))
				s = s.concat(c + "");
		}
		return s;
	}

}
