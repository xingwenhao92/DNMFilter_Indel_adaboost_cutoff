package edu.duke.humangenome.file;

public class VisualData {
	private String editcigar;
	private String direction;
	private int start;
	private int end;
	private int numb;//序号
	private String cigar;
	//private String read;
	private String individualID;
	private int MappingQuality;
	private String ReadBases;
	public VisualData(boolean direction, String cigar, int start, int end, String readbases, int MappingQuality, int numb, String individualID) {
		this.cigar = cigar;//edit cigar for html easy to handle change to format R:1;H:1;M:50;S:40 H:1;M:50;S:40;F:1  F表示正链放末尾 R表示反链放开头
		this.direction = (direction == false) ? "positive" : "negative";
		this.start = start;
		this.end = end;
		this.ReadBases = readbases;
		this.MappingQuality = MappingQuality;
		this.numb = numb;
		this.individualID = individualID;
		this.editcigar = formatCigar(cigar.trim(), this.direction);
	}
	
	public String formatCigar(String cigar, String direction) {
		String [] notnumber = cigar.split("[0-9]+");
		String [] number = cigar.split("[^0-9]+");
		StringBuilder strb = new StringBuilder();
		//System.out.println(number.length + "\t" + notnumber.length);
		/*for (int i = 0; i < number.length; i++) {
			System.out.println(number[i]);
		}
		//System.out.println();
		for(int i = 0; i < notnumber.length; i++) {
			System.out.println(notnumber[i]);
		}
		//System.out.println("\n");*/
		if (direction.equals("negative")) {
			strb.append("R:1;");
		}
		
		for (int i = 0; i < number.length - 1; i++) {
			strb.append(notnumber[i + 1] + ":" + number[i] + ";");
		}
		strb.append(notnumber[number.length] + ":" + number[number.length - 1]);
		if (direction.equals("positive")) {
			strb.append(";F:1");
		}
		return strb.toString();
	}
	
	public String getEditCigar() {
		return editcigar;
	}
	
	public String getDirection() {
		return direction;
	}
	
	public String getIndividualID() {
		return individualID;
	}
	
	public int getStart() {
		return start;
	}
	
	public int getNumb() {
		return numb;
	}
	
	public int getEnd() {
		return end;
	}
	
	public String getCigar() {
		return cigar;
	}
	
	public int getMappingQuality() {
		return MappingQuality;
	}
	
	public String getRead() {
		return ReadBases;
	}
}
