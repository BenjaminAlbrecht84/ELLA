package ella.model.aligner.utils;

public enum MyStrand {

	FWD("+"), REV("-");

	private final String id;

	private MyStrand(String s) {
		id = s;
	}

	public String getID() {
		return id;
	}

	public int getNumericalID() {
		return id.equals("+") ? +1 : -1;
	}

}
