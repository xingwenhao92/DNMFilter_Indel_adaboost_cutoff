/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.test;

import edu.duke.humangenome.sam.IndelPileup;
import java.io.File;
import java.util.Iterator;
import java.util.List;
import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.reference.ReferenceSequenceFileWalker;
import net.sf.picard.util.Interval;
import net.sf.picard.util.IntervalList;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.util.SamLocusIterator.RecordAndOffset;
import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;

/**
 *
 * @author Yongzhuang Liu
 */
public class Test {

    public void test() {

        String filename = "D:/Data/Chr1.Test.bam";
        
        String referenceSequenceFile="D:/Data/hs37d5.fa";
        ReferenceSequenceFileWalker referenceSequenceFileWalker = new ReferenceSequenceFileWalker(new File(referenceSequenceFile));

        ReferenceSequence referenceSequence = referenceSequenceFileWalker.get(0);

        SAMFileReader samFileReader = new SAMFileReader(new File(filename));
        samFileReader.setValidationStringency(SAMFileReader.ValidationStringency.SILENT);

        IntervalList intervalList = new IntervalList(samFileReader.getFileHeader());
        intervalList.add(new Interval("1", 61871, 61871));
        SamLocusIterator samLocusIterator = new SamLocusIterator(samFileReader, intervalList);

        Iterator<SamLocusIterator.LocusInfo> iterator = samLocusIterator.iterator();
        while (iterator.hasNext()) {
            SamLocusIterator.LocusInfo locusInfo = iterator.next();
            String refAllele = "C";
            String altAllele = "CT";
            IndelPileup indelPileup=new IndelPileup(referenceSequence, locusInfo, refAllele, altAllele);
            System.out.println("Ref Allele Count: "+ indelPileup.getRefAlleleCount());
            System.out.println("Alt Allele Count: "+ indelPileup.getAltAlleleCount());
            
        }
    }


    public static void main(String[] args){  
        (new Test()).test();    
//         System.out.println( CigarOperator.I.consumesReadBases());
//         System.out.println( CigarOperator.I.consumesReferenceBases());
        
        
    }
}
