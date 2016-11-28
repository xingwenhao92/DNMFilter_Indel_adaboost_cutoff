/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package edu.duke.humangenome.core;

import edu.duke.humangenome.sam.MultiPileup;
import edu.duke.humangenome.sam.IndelPileup;
import edu.duke.humangenome.util.BaseUtils;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Map;
import java.util.Properties;
import java.util.regex.Matcher;
import java.util.regex.Pattern;

import net.sf.picard.reference.ReferenceSequence;
import net.sf.picard.util.SamLocusIterator;
import net.sf.picard.reference.IndexedFastaSequenceFile;
/**
 *
 * @author Yongzhuang Liu
 */
public class FeatureSelection {

    private final static int FATHER_INDEX = 0;
    private final static int MOTHER_INDEX = 1;
    private final static int OFFSPRING_INDEX = 2;
    private ReferenceSequence referenceSequence;
    private MultiPileup trioPileup;
    private Properties properties;
    private IndexedFastaSequenceFile indexfasta;
    public FeatureSelection(ReferenceSequence referenceSequence, MultiPileup trioPileup, IndexedFastaSequenceFile indexfasta) {
        this.referenceSequence = referenceSequence;
        this.trioPileup = trioPileup;
        this.indexfasta = indexfasta;
        
    }

    public String extract() throws FileNotFoundException {
    	//System.out.println(properties.getProperty("reference"));
    	//indexfasta = new IndexedFastaSequenceFile(new File(properties.getProperty("reference")));
        String chrom = trioPileup.getReferenceName();
        //System.out.println("chrom"+chrom);
        int pos = trioPileup.getPosition();
        int homostartpos = pos - 5;
        int homoendpos = pos + 5;
        int strstartpos = pos - 15;
        int strendpos = pos + 15;
        //int midpos = (strendpos - strstartpos)/2;
        //String homo = new String(indexfasta.getSubsequenceAt(chrom, homostartpos, homoendpos).getBases()); 
        String str = new String(indexfasta.getSubsequenceAt(chrom, strstartpos, strendpos).getBases());
        String homo = str.substring(homostartpos-strstartpos, homoendpos-strstartpos+1);
        //System.out.println("pos"+pos);
//        char ref = getRefBase();
//        char alt = getMostAltAllele();
        int[] index = new int[]{FATHER_INDEX, MOTHER_INDEX, OFFSPRING_INDEX};

        StringBuilder recordBuilder = new StringBuilder();

        int[] fAlleleCount = new int[2];
        int[] mAlleleCount = new int[2];
        int[] oAlleleCount = new int[2];

        for (int i = 0; i < index.length; i++) {
            
            IndelPileup pileup = new IndelPileup(referenceSequence,trioPileup.getLocusInfo(index[i]),trioPileup.getRefAllele(), trioPileup.getAltAllele());
            int readDepth = pileup.getDepth();

//            int meanRefBaseQuality = pileup.getMeanRefMappingQuality();
//            int meanAltBaseQuality = pileup.getMeanAltMappingQuality();

            int meanRefMappingQuality = pileup.getMeanRefMappingQuality();
            int meanAltMappingQuality = pileup.getMeanAltMappingQuality();

            int meanRefDistanceToThreePrime = pileup.getMeanRefDistanceToThreePrime();
            int meanAltDistanceToThreePrime = pileup.getMeanAltDistanceToThreePrime();

            double refFractionOfMQ0Reads = pileup.getRefFractionOfMQ0Reads();
            double altFractionOfMQ0Reads = pileup.getAltFractionOfMQ0Reads();

            double refFractionOfSoftClippedReads = pileup.getRefFractionOfSoftClippedReads();
            double altFractionOfSoftClippedReads = pileup.getAltFractionOfSoftClippedReads();

            int forwardRef = pileup.getForwardRefCount();
            int reverseRef = pileup.getReverseRefCount();
            int forwardAlt = pileup.getForwardAltCount();
            int reverseAlt = pileup.getReverseAltCount();

            int strandBias = pileup.getStrandBias();

            int refStrandDirection;
            int altStrandDirection;
            if (pileup.getRefStrandDirection()) {
                refStrandDirection = 1;
            } else {
                refStrandDirection = 0;
            }
            if (pileup.getAltStrandDirection()) {
                altStrandDirection = 1;
            } else {
                altStrandDirection = 0;
            }

            if (index[i] == FATHER_INDEX) {
                fAlleleCount[0] = forwardRef + reverseRef;
                fAlleleCount[1] = forwardAlt + reverseAlt;
            }
            if (index[i] == MOTHER_INDEX) {
                mAlleleCount[0] = forwardRef + reverseRef;
                mAlleleCount[1] = forwardAlt + reverseAlt;
            }
            if (index[i] == OFFSPRING_INDEX) {
                oAlleleCount[0] = forwardRef + reverseRef;
                oAlleleCount[1] = forwardAlt + reverseAlt;
            }
            
            double alleleBalance = 0;//防止出现除零错误
            if( (pileup.getRefAlleleCount() + pileup.getAltAlleleCount()) != 0 )
            	alleleBalance = (double) pileup.getAltAlleleCount() / (pileup.getRefAlleleCount() + pileup.getAltAlleleCount());
            	
            /*
             double meanRefNearbyIndels = (double) pileup.getNearbyIndels(ref) / (forwardRef + reverseRef);
             double meanAltNearbyIndels = (double) pileup.getNearbyIndels(alt) /  (forwardAlt+ reverseAlt);
             double meanRefNearbyMismatches = (double) (pileup.getOverallMismatches(referenceSequence, ref)) / (forwardRef + reverseRef);
             double meanAltNearbyMismatches = (double) (pileup.getOverallMismatches(referenceSequence, alt) - forwardAlt - reverseAlt) / (forwardAlt+ reverseAlt);

             * 
             */


            double meanRefNearbyIndels = pileup.getMeanRefNearbyIndels();
            double meanAltNearbyIndels = pileup.getMeanAltNearbyIndels();
            double meanRefNearbyMismatches = pileup.getMeanRefNearbyMismatches();
            double meanAltNearbyMismatches = pileup.getMeanAltNearbyMismatches();

            recordBuilder.append(alleleBalance + ","  + readDepth + ",");
            recordBuilder.append(meanRefMappingQuality + "," + meanAltMappingQuality + ",");
            recordBuilder.append(meanRefDistanceToThreePrime + "," + meanAltDistanceToThreePrime + ",");
            recordBuilder.append(refFractionOfMQ0Reads + "," + altFractionOfMQ0Reads + ",");
            recordBuilder.append(refFractionOfSoftClippedReads + "," + altFractionOfSoftClippedReads + ",");
            recordBuilder.append(meanRefNearbyMismatches + "," + meanAltNearbyMismatches + ",");
            recordBuilder.append(meanRefNearbyIndels + "," + meanAltNearbyIndels + ",");
            recordBuilder.append(refStrandDirection + "," + altStrandDirection + ",");
            recordBuilder.append(strandBias + ",");
        }
        //System.out.println(homo + "\t" + str);
        //System.out.println(gethomoflag(homo) + "\t" + getstrflag(str));
        int homopolymerflag = gethomoflag(homo);
        int shorttandemrepeatflag = getstrflag(str);
        int pValueOfFatherToOffspring = getPhredPValue(fAlleleCount, oAlleleCount);
        int pValueOfMotherToOffspring = getPhredPValue(mAlleleCount, oAlleleCount);
        recordBuilder.append(pValueOfFatherToOffspring + "," + pValueOfMotherToOffspring + "," + homopolymerflag + "," +shorttandemrepeatflag);///加入了homo和short tandem repeat
        
        return recordBuilder.toString();
    }
    
    
    public int getPhredPValue(int[] pAlleleCount, int[] oAlleleCount) {
        FisherExact fe = new FisherExact(pAlleleCount[0] + pAlleleCount[1] + oAlleleCount[0] + oAlleleCount[1]);
        double p = fe.getCumlativeP(pAlleleCount[0], pAlleleCount[1], oAlleleCount[0], oAlleleCount[1]);
        if (p > 1) {
            p = 1.0;
        }
        return (int) (Math.log10(p) * (-10));
    }
    
    public int gethomoflag(String ahomo) {
    	Pattern pa;
    	Matcher ma;
    	pa = Pattern.compile("([atcg])\\1{4,}",Pattern.CASE_INSENSITIVE);
    	ma = pa.matcher(ahomo);
    	while (ma.find()) 
    		return 1;
    	return 0;
    }
    
    public int getstrflag(String astr) {
    	int[][] motifandtimes = {{2,3},{3,3},{4,3},{5,3}};
    	Pattern pa;
    	Matcher ma;
    	for (int i = 0; i <motifandtimes.length; i++) {
    		pa = Pattern.compile("([atcg]{" + motifandtimes[i][0] + "})\\1{" + (motifandtimes[i][1] - 1) + ",}",Pattern.CASE_INSENSITIVE);
			ma = pa.matcher(astr);
			while (ma.find()) {
				String substr = ma.group();
				int start = ma.start();
				int end = ma.end();//start 是索引的值从零开始包含start  end 是偏移从1开始…… 包含end
				if (start <= (astr.length()-1)/2+1 && end >= (astr.length()+1)/2-1) {
					//System.out.println(start+ " " + end);
					if (isHomo(ma.group().toLowerCase(), motifandtimes[i][0]) == false) {
						//writer.println(titleline + "\t" + strline + "\t" +start + "\t" + end);
						//writer.println(strline);
						return 1;
					}
				}
			}
    	}
    	return 0;
    }
    
    boolean isHomo(String str, int motiflen) {
		//System.out.println(str);
		char achr = 0;
		if (str.length() != 0) 
			achr = str.charAt(0);
		for (int i = 1; i < motiflen && i < str.length(); i++) {
			if (str.charAt(i) != achr)
				return false;
		}
		return true;
	}
}
