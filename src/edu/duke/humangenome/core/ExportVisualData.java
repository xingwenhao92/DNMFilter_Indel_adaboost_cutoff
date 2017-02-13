package edu.duke.humangenome.core;

import java.io.File;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Properties;

import org.apache.log4j.Logger;

import edu.duke.humangenome.file.ConfigReader;
import edu.duke.humangenome.file.DNMReader;
import edu.duke.humangenome.file.DNMRecord;
import edu.duke.humangenome.file.PEDReader;
import edu.duke.humangenome.file.Trio;
import edu.duke.humangenome.file.VisualData;
import net.sf.picard.reference.IndexedFastaSequenceFile;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

public class ExportVisualData {
	private static Logger logger = Logger.getLogger(ExportVisualData.class);
    private final static int FATHER_INDEX = 0;
    private final static int MOTHER_INDEX = 1;
    private final static int OFFSPRING_INDEX = 2;
	private Properties properties;
	String outputPath;
	Properties bams;
	File referenceSequenceFile;
	List<Trio> trios;
	String positon;
	private IndexedFastaSequenceFile indexfasta;
	public ExportVisualData(Properties properties) throws IOException {
		this.properties = properties;
		this.referenceSequenceFile = new File(properties.getProperty("reference"));
        this.trios = (new PEDReader(properties.getProperty("pedigree"))).getTrios();
        this.bams = (new ConfigReader(properties.getProperty("bam"))).parse();
        this.outputPath = new File(properties.getProperty("outputpath")).getAbsolutePath();
        this.positon = properties.getProperty("position");
        this.indexfasta = new IndexedFastaSequenceFile(new File(properties.getProperty("reference")));
	}
	
	public void run() throws IOException {//seperate by position but not trios so this will not use the structure in GBM.java
		DNMReader dnmposition = new DNMReader(properties.getProperty("position"));
		
		HashMap<String, Trio> map = new HashMap();
		for(int i = 0; i < trios.size(); i++) {
			map.put(trios.get(i).getFamilyID(), trios.get(i));
		}
		
		DNMRecord dnmRecord = null;
		dnmRecord = dnmposition.getNextRecord();
		String familyID = null;
		
		while (dnmRecord != null) {
			//System.out.println(outputPath);
			String folder = outputPath + "/";
			familyID = dnmRecord.getFamilyID();
			Trio trio = map.get(familyID);
			String[] individualID = new String[3];
			individualID[0] = trio.getFather().getIndividualID(); 
			individualID[1] = trio.getMother().getIndividualID();
			individualID[2] = trio.getOffspring().getIndividualID();
									
			SAMFileReader[] trioSAMFileReader = new SAMFileReader[3];
            //IntervalList[] trioIntervalList = new IntervalList[3];
            SamLocusIterator[] samLocusIterator = new SamLocusIterator[3];
            trioSAMFileReader[0] = new SAMFileReader(new File(bams.getProperty(trio.getFather().getIndividualID())));
            trioSAMFileReader[1] = new SAMFileReader(new File(bams.getProperty(trio.getMother().getIndividualID())));
            trioSAMFileReader[2] = new SAMFileReader(new File(bams.getProperty(trio.getOffspring().getIndividualID())));
            for (int i = 0; i < 3; i++) {
                trioSAMFileReader[i].setValidationStringency(SAMFileReader.ValidationStringency.SILENT);
                //trioIntervalList[i] = new IntervalList(trioSAMFileReader[i].getFileHeader());
            }
            
            //Interval interval = new Interval(dnmRecord.getChrom(), dnmRecord.getPos(),dnmRecord.getPos());
            List<VisualData> visualdataList = new ArrayList<>();
            SAMRecordIterator iter = null;
            int minstart = Integer.MAX_VALUE;
            int maxend = Integer.MIN_VALUE;
            for (int i = 0; i < 3; i++) {
            	iter = trioSAMFileReader[i].queryOverlapping(dnmRecord.getChrom(), dnmRecord.getPos(), dnmRecord.getPos());
            	int count = 0;
            	
            	while (iter.hasNext()) {
            		count++;
            		SAMRecord SamRead = iter.next();
            		String cigar = SamRead.getCigarString();
            		int start  = SamRead.getAlignmentStart();
            		int end = SamRead.getAlignmentEnd();
            		String readbases = SamRead.getReadString();
            		int mappingquality = SamRead.getMappingQuality();
            		boolean direction = SamRead.getReadNegativeStrandFlag();
            		//System.out.println(start + " " + end);
            		visualdataList.add(new VisualData(direction, cigar, start, end, readbases, mappingquality, count, individualID[i]));
            		//if (start < end) {
            		minstart = Integer.min(minstart, start);
            		maxend = Integer.max(maxend, end);
            		//} else {
            		//	minstart = Integer.min(minstart, end);
            		//	maxend = Integer.max(maxend, end);
            		//}
            	}
            	iter.close();
            }
            minstart--;
            maxend++;
            String refread = new String(indexfasta.getSubsequenceAt(dnmRecord.getChrom(), minstart, maxend).getBases());
            File tempfile = new File(folder + familyID + "_" + dnmRecord.getChrom() + "_" + dnmRecord.getPos() + ".csv");
            File reffile = new File(folder + familyID + "_" + dnmRecord.getChrom() + "_" + dnmRecord.getPos() + ".csv.ref");
            PrintWriter writer = new PrintWriter(tempfile);
            PrintWriter refwriter = new PrintWriter(reffile);
            //writer.println("#reference," + minstart + "," + maxend + "," + refread);
            //writer.println("#position," + dnmRecord.getChrom() + "," + dnmRecord.getPos());
            VisualData tmp = null;
            writer.println("allID,groupID,groupName,Direction,startPosition,endPosition,mappingquality,cigar,editcigar,baseread");
            int countall = 0;
            for (int i = 0; i < visualdataList.size(); i++) {
            	tmp = visualdataList.get(i);
            	countall ++;
            	/*if (tmp.getNumb() == 1) {
            		writer.println("#" + tmp.getIndividualID());
            	}*/
            	writer.println(countall + "," + tmp.getNumb() + "," + tmp.getIndividualID() +  "," + tmp.getDirection() + "," + tmp.getStart() + "," + tmp.getEnd() + ","+ tmp.getMappingQuality() + "," + tmp.getCigar() + "," + tmp.getEditCigar() + "," + tmp.getRead());
            }
            refwriter.println("chrom,visualPositon,startPosition,endPosition,refread");
            refwriter.println(dnmRecord.getChrom() + "," + dnmRecord.getPos() + "," + minstart + "," + maxend + "," + refread);
            writer.close();
            refwriter.close();
            dnmRecord = dnmposition.getNextRecord();
		}		
	}
}
