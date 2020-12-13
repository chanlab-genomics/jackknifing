import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;



public class Jackknife {
	
	/**
	 * List that contains sequences.
	 */
	static ArrayList<String> listSeq = new ArrayList<String>();
	
	/**
	 * List that contains names of the sequences.
	 */
	static ArrayList<String> seqsNames = new ArrayList<String>();
	
	/**
	 * Function used to print the DistanceMatrix in PHYLIP format in a txt file.
	 * @param matCor	normalized matrix 
	 * @param name	List of sequences name
	 * @param path	Path of the fasta file (input)
	 * @param k	parameter k used
	 * 
	 * @since 1.0
	 */
	public static void printMatrixPhylip(ArrayList<ArrayList<String>> pseudomat, String path){
		
		PrintWriter f = null;
		//String[] new_path = path.split(".fas");
		//new_path[0]+"_k"+k+"_"+d2fam+".txt"
		try {
			f  = new PrintWriter(new BufferedWriter(new FileWriter(path)));
		} catch (IOException e) {
			e.printStackTrace();
		}
		f.print("     "+pseudomat.size()+"\n");
		
		int name_sizeOne=pseudomat.get(0).get(0).length();
		String spaceOne ="";
		if (name_sizeOne<10){
			for (int tmp=name_sizeOne;tmp<10;tmp++){
				spaceOne+=" ";
			}
			pseudomat.get(0).set(0, pseudomat.get(0).get(0)+spaceOne);
		}
		if (pseudomat.get(0).get(0).length() > 47){
			pseudomat.get(0).set(0, pseudomat.get(0).get(0).substring(0, 46)+"...");
		}
		
		f.print(pseudomat.get(0).get(0)+"      \n");
		
		int size_name=pseudomat.size();
		for (int i=1; i<size_name;i++){
			
			int name_size=pseudomat.get(i).get(0).length();
			if (name_size>=10){
				pseudomat.get(i).set(0, pseudomat.get(i).get(0).substring(0, 9)+" ");
			}
			else{
				String space ="";
				for (int tmp=name_size;tmp<10;tmp++){
					space+=" ";
				}
				pseudomat.get(i).set(0, pseudomat.get(i).get(0)+space);
			}
			
			
			f.print(pseudomat.get(i).get(0));
			for (int j=1; j<pseudomat.get(i).size();j++){
				
				String tmp = pseudomat.get(i).get(j);
				f.print("  "+tmp);
			}
			f.print(" \n");
		}
		
		f.close();
	}
	
	/**
	 * Function used to:
	 * 1: generate the readerFasSeq object (call readerTest)
	 * 2: complete the listSeq and seqsNames
	 * Steps 2 is done by the readerTest object (see ReaderFasSeq.java)
	 * @param path	fasta or fastq file path.
	 * 
	 * @see ReaderFasSeq#readFromFol()
	 * 
	 * @since 1.0
	 */
	public static void readFolder(String path){
		File di   = new File(path);
		File fl[] = di.listFiles();
		Loader readerTest;
		for (int i=0; i < fl.length; i++)
		{
			
			System.out.println("filename: "+fl[i].toString());
			readerTest = new Loader(fl[i].toString());
			if (readerTest.loadFasta(fl[i].toString())){
				ArrayList<String>  list1 = new ArrayList<String>();
				list1 = readerTest.readFromFol();
				listSeq.add(list1.get(0));
				seqsNames.add(list1.get(1));
			}
			else{}
		}
	}
	
	public static void write_seqs(String path, ArrayList<String> m_seqs){
		
		File di   = new File(path);
		File fl[] = di.listFiles();
		PrintWriter f =null;
		Loader readerTest;
		int count=0;
		
		for (int i=0; i < fl.length; i++)
		{
			readerTest = new Loader(fl[i].toString());
			if (readerTest.loadFasta(fl[i].toString())){
				try {
					String tmp_file = fl[i].toString();
					String[] new_path = tmp_file.split(".fna");
					f = new PrintWriter(new BufferedWriter(new FileWriter(new_path[0]+"_rep.fna")));
					f.write(">"+seqsNames.get(count));
					f.write("\n");
					f.write(m_seqs.get(count));
					f.write("\n");
					f.write("\n");
					count++;
				} catch (IOException e) {
					// TODO Auto-generated catch block
					e.printStackTrace();
				}
				
				f.close();
			}
		}
	}
	
	public static void main(String[] args) {
		// TODO Auto-generated method stub
		String path=args[0];
		readFolder(path);
		JackKnifing Jackknifer;
		Jackknifer = new JackKnifing(listSeq);
		ArrayList<String> pseudoSeqs = new ArrayList();
		int portion = 100 - Integer.parseInt(args[1]);
		pseudoSeqs = Jackknifer.datasampling(portion);
		write_seqs(path,pseudoSeqs);
		
	}

}
