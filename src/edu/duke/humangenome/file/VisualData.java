package edu.duke.humangenome.file;

import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.samtools.AlignmentBlock;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMRecord;

public class VisualData {
	private String editcigar;
	private String direction;
	private int start;
	private int startContainSoftAndHardClip;
	private int endContainSoftAndHardClip;
	private int end;
	private int numb;//序号
	private String mismatch;
	private String cigar;
	//private String read;
	private String individualID;
	private int MappingQuality;//if mapping quality is zero then set blank
	private String ReadBases;
	public VisualData(DNMRecord dnmRecord, IndexedFastaSequenceFile indexfasta, SAMRecord SamRead, int numb, String individualID) {
		this.cigar = SamRead.getCigarString();//edit cigar for html easy to handle change to format R:1;H:1;M:50;S:40 H:1;M:50;S:40;F:1  F表示正链放末尾 R表示反链放开头
		this.direction = (SamRead.getReadNegativeStrandFlag() == false) ? "positive" : "negative";
		this.start = SamRead.getAlignmentStart();
		this.end = SamRead.getAlignmentEnd();
		this.startContainSoftAndHardClip = start;
		this.endContainSoftAndHardClip = end;
		this.ReadBases = SamRead.getReadString();
		this.MappingQuality = SamRead.getMappingQuality();
		this.numb = numb;
		this.individualID = individualID;
		this.editcigar = formatCigar(SamRead, this.direction);
		this.mismatch = findMismatch(SamRead, indexfasta, dnmRecord);
	}
	
	public String formatCigar(SAMRecord SamRead, String direction) {
		StringBuilder strb = new StringBuilder();
		int k = 0;
		if (direction.equals("negative")) {
			strb.append("R:1;");
		}
		
		for (CigarElement e : SamRead.getCigar().getCigarElements()) {
			k ++;
			String operat = e.getOperator().toString();
			int operatlen = e.getLength();
			if (operat.equals("S")) { // operat.equals("H") do not add to length   k = 1 "S" in the first place  k != 1 "S" in the last place
				if (k == 1)startContainSoftAndHardClip = start - operatlen;
				else endContainSoftAndHardClip = end + operatlen;
			}
			strb.append(operat + ":" + operatlen + ";");
		}
		if (direction.equals("positive")) {
			strb.append("F:1;");
		}
		
		/*String [] notnumber = cigar.split("[0-9]+");
		String [] number = cigar.split("[^0-9]+");		
		//System.out.println(number.length + "\t" + notnumber.length);
		for (int i = 0; i < number.length; i++) {
			System.out.println(number[i]);
		}
		//System.out.println();
		for(int i = 0; i < notnumber.length; i++) {
			System.out.println(notnumber[i]);
		}
		//System.out.println("\n");
		
		//calculate end and start 
		if (notnumber[1].equals("H") || notnumber[1].equals("S")) {
			startContainSoftAndHardClip = start - Integer.parseInt(number[0]);
		}
		if(notnumber[number.length].equals("H") || notnumber[number.length].equals("S")) {
			endContainSoftAndHardClip = end + Integer.parseInt(number[number.length - 1]);
		}
		
		if (direction.equals("negative")) {
			strb.append("R:1;");
		}
		
		//if (MappingQuality == 0) //mapq = 0 then set to white
		
		for (int i = 0; i < number.length - 1; i++) {
			strb.append(notnumber[i + 1] + ":" + number[i] + ";");
		}
		strb.append(notnumber[number.length] + ":" + number[number.length - 1]);

		if (direction.equals("positive")) {
			strb.append(";F:1");
		}*/
		return strb.toString().substring(0, strb.length()-1);//remove the last ";"
	}
	
	public String findMismatch(SAMRecord SamRead, IndexedFastaSequenceFile indexfasta, DNMRecord dnmRecord) {
		StringBuilder misstr = new StringBuilder();// the foramt is A:1232;T:1235
		String baseRead = SamRead.getReadString();
		for (AlignmentBlock e : SamRead.getAlignmentBlocks()) {
			String refRead = new String(indexfasta.getSubsequenceAt(dnmRecord.getChrom(), e.getReferenceStart(), e.getReferenceStart()+e.getLength()-1).getBases());
			for (int j = 0; j < e.getLength(); j++) {
				if (baseRead.charAt(e.getReadStart() + j - 1) != refRead.charAt(j)) {
					misstr.append(baseRead.charAt(e.getReadStart() + j - 1) + ":" + (e.getReferenceStart() + j) + ";");
				}
			}
			//System.out.println(e.getLength() + "\t" + e.getReadStart() + "\t" + e.getReferenceStart() + "\t" + refRead + "\t" + baseRead.substring(e.getReadStart()-1 , e.getReadStart()+e.getLength()-1));
		}
		/*if (misstr.length() != 0)
			System.out.println(misstr.toString().substring(0, misstr.length()-1));*/
		if (misstr.length() == 0) {
			return "none";
		}
		return misstr.toString().substring(0, misstr.length()-1);		
	}
	
	public String getEditCigar() {
		return editcigar;
	}
	public String getMismatch() {
		return mismatch;
	}
	
	public String getDirection() {
		return direction;
	}
	
	public String getIndividualID() {
		return individualID;
	}
	
	public int getstartContainSoftAndHardClip() {
		return startContainSoftAndHardClip;
	}
	
	public int getendContainSoftAndHardClip() {
		return endContainSoftAndHardClip;
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
