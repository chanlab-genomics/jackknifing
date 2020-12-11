import java.io.BufferedReader;
import java.io.FileNotFoundException;
import java.io.FileReader;
import java.io.IOException;
import java.io.Reader;
import java.util.ArrayList;


public class Loader {

	private Reader doc;//file Reader.
	private String file;
	private String alphabet = "acgt";
	
	public Loader(String m_path){
		this.file = m_path;
	}//constructor
	
	/**
	 * Function used to load a FASTA file and verified if the path is a good path
	 * @param fileName
	 * @return Boolean	boolean to verify the path.
	 * 
	 * @since 1.0
	 */
	public boolean loadFasta(String fileName){
		
		try{
			if (fileName.endsWith(".ffn") || fileName.endsWith(".fna") || fileName.endsWith(".faa") || fileName.endsWith(".frn") || fileName.endsWith(".fas") || fileName.endsWith(".fasta")){
				this.doc = new FileReader(fileName);
				return true;
			}
			else{
				return false;
			}
		}
		
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	        return false;
		}
	        
	    catch (@SuppressWarnings("hiding") IOException e) {
	        e.printStackTrace();
	        return false;
	    }
		
	}
	
	/**
	 * Function used to read a FASTA file
	 * 1: stock sequences in listSeq List.
	 * 2: stock sequences name in name List.
	 * 3: return these 2 List in the Main.
	 * @return Seqs_names	list of sequence.
	 * 
	 * @since 1.0
	 */
	public ArrayList<ArrayList<String>> readFasDoc(){
		
		ArrayList<String> listSeq = new ArrayList<String>();
		ArrayList<String> name = new ArrayList<String>();
		BufferedReader loadBuffer = null;
		String line;
		String sequence = new String();
		loadBuffer = new BufferedReader(this.doc);
		
		try {
			while ((line = loadBuffer.readLine()) != null){
				try{
					if (line.charAt(0) == '>'){
						String tmp[] = line.split(">");
						String name2 = tmp[1];
						name.add(name2);
						if (sequence.length()>0){
							listSeq.add(sequence);
							sequence = "";
						}
					}
					
					if (line.charAt(0) != '>' && line.charAt(0) != '!'){
						sequence += line;
					}
				}
				catch(StringIndexOutOfBoundsException e){}
			}
			listSeq.add(sequence);
			loadBuffer.close();
		} 
		
		catch (IOException e) {
			e.printStackTrace();
		}
		
		ArrayList<ArrayList<String>> Seqs_names = new ArrayList<ArrayList<String>>();
		Seqs_names.add(listSeq);
		Seqs_names.add(name);
		return Seqs_names;
		
	}
	
	/**
	 * Function used to load a matrix file and verified if the path is a good path
	 * @param fileName
	 * @return Boolean	boolean to verify the path.
	 * 
	 * @since 1.0
	 */
	public boolean loadMatrix(String fileName){
		
		try{
			if (fileName.endsWith(".txt")){
				this.doc = new FileReader(fileName);
				return true;
			}
			else{
				return false;
			}
		}
		
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	        return false;
		}
	        
	    catch (@SuppressWarnings("hiding") IOException e) {
	        e.printStackTrace();
	        return false;
	    }
		
	}
	
	/**
	 * Function used to read a FASTA file from a folder
	 * 1: stock sequences in listSeq List.
	 * 2: stock sequences name in name List.
	 * 3: return these 2 List in the Main.
	 * @return Seqs_names	list of sequence.
	 * 
	 * @since 1.0
	 */
	public ArrayList<String> readFromFol(){
		/*
		ArrayList<String> Seqs_names = null;
		if (alphabet.equals("acgt")){
			LinkedHashMap<String, DNASequence> tmp=null;
			try {
				tmp = FastaReaderHelper.readFastaDNASequence(new File(file));
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
				
			}
			
			Seqs_names = new ArrayList<String>();
			for (Entry<String, DNASequence> tmp3 : tmp.entrySet()){
				
				Seqs_names.add(tmp3.getValue().getSequenceAsString());
				Seqs_names.add(tmp3.getKey());
				break;
			}
		}
		
		else{
			LinkedHashMap<String, ProteinSequence> tmp2=null;
			try {
				tmp2 = FastaReaderHelper.readFastaProteinSequence(new File(file));
			} catch (Exception e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}
			
			Seqs_names = new ArrayList<String>();
			for (Entry<String, ProteinSequence> tmp3 : tmp2.entrySet()){
				
				Seqs_names.add(tmp3.getValue().getSequenceAsString());
				Seqs_names.add(tmp3.getKey());
				break;
			}
		}*/
		
		
		String sequence = new String();
		BufferedReader loadBuffer = null;
		String line;
		String name = null;
		loadBuffer = new BufferedReader(this.doc);
		
		try {
			line = loadBuffer.readLine();
			if (line.charAt(0) == '>'){
				String tmp[] = line.split(">");
				name = tmp[1];
			}
			System.out.println(name);
			while ((line = loadBuffer.readLine()) != null){
				try{
					if (line.charAt(0) != '>' && line.charAt(0) != '!'){
						sequence += line;
					}
					else{
						//System.out.println(line);
					}
				}
				catch(StringIndexOutOfBoundsException e){}
			}
			loadBuffer.close();
		} 
		
		catch (IOException e) {
			e.printStackTrace();
		}
		ArrayList<String> Seqs_names = new ArrayList<String>();
		//System.out.println(sequence);
		
		Seqs_names.add(sequence);
		Seqs_names.add(name);
		
		return Seqs_names;
		
	}
	
	/**
	 * Function used to convert sequences in lower case and delete the incorrect characters.
	 * @param listSeq	list of sequences.
	 * @return n_listSeq	list of sequences.
	 * 
	 * @since 1.0
	 */
	public ArrayList<String> convertSequences(ArrayList<String> listSeq){
		
		ArrayList<String> n_listSeq = new ArrayList<String>();
		int N = listSeq.size();
		for (int i = 0; i<N;i++){
			String Seq_tmp=listSeq.get(i).replaceAll("\\*","");
			if (alphabet.equals("acgt")){
				Seq_tmp=Seq_tmp.replaceAll("[U,W,S,M,K,R,Y,B,D,H,V,N,X]","");
			}
			
			else if (alphabet.equals("acgu")){
				Seq_tmp=Seq_tmp.replaceAll("[T,W,S,M,K,R,Y,B,D,H,V,N,X]","");
			}
			
			else {
				Seq_tmp=Seq_tmp.replaceAll("[B,X]","");
			}
			n_listSeq.add(Seq_tmp.toLowerCase());
		}
		
		return n_listSeq;
	}
	
	
	public boolean loadTree(String fileName){
		
		try{
			if (fileName.endsWith(".tre")){
				this.doc = new FileReader(fileName);
				return true;
			}
			else{
				return false;
			}
		}
		
		catch (FileNotFoundException e) {
	        e.printStackTrace();
	        return false;
		}
	        
	    catch (@SuppressWarnings("hiding") IOException e) {
	        e.printStackTrace();
	        return false;
	    }
		
	}
	
	/**
	 * Function used to read a matrix file
	 * 1: stock sequences in listSeq List.
	 * 2: stock sequences name in name List.
	 * 3: return these 2 List in the Main.
	 * @return Seqs_names	list of sequence.
	 * 
	 * @since 1.0
	 */
	public ArrayList<ArrayList<String>> readMatrix(){
		
		ArrayList<ArrayList<String>>  matrix = new ArrayList<ArrayList<String>> ();
		BufferedReader loadBuffer = null;
		String line;
		loadBuffer = new BufferedReader(this.doc);
		
		int number;
		try {
			//number of species
			line = loadBuffer.readLine();
			line = line.replaceAll(" ", "");
			number = Integer.parseInt(line);
			//System.out.println(number);
			
			//first species (only name)
			line = loadBuffer.readLine();
			line = line.replaceAll("   ",";");
			line = line.replaceAll("  ",";");
			String[] tmp = line.split(";");
			//System.out.println(tmp[0]);
			ArrayList<String> row1 = new ArrayList<String> ();
			row1.add(tmp[0]);
			matrix.add(row1);
			
			while ((line = loadBuffer.readLine()) != null){
				try{
					line = line.replaceAll("   ",";");
					line = line.replaceAll("  ",";");
					String[] tmp2 = line.split(";");
					//System.out.println(tmp2[1]);
					ArrayList<String> rowx = new ArrayList<String> ();
					for (int i=0; i<tmp2.length;i++){
						rowx.add(tmp2[i]);
					}
					matrix.add(rowx);
				}
				
				catch(StringIndexOutOfBoundsException e){}
			}
		}
		
		catch (IOException e) {
			e.printStackTrace(); 
		}
		
		System.out.println(matrix);
		return matrix;
		
	}
	
	public ArrayList<ArrayList> readTree(){
		
		Tree tree1 = new Tree();
		BufferedReader loadBuffer = null;
		String line;
		int n_node=0;
		int n_leaf=0;
		int father_node=0;
		String label = "";
		String b_length = "";
		int check=0;
		
		loadBuffer = new BufferedReader(this.doc);
		try {
			line = loadBuffer.readLine();
			char[] line2 = line.toCharArray();
			int N = line2.length;
			for (int i=0;i<N;i++){
				if (line2[i]=='('){
						n_node++;
						Node tmp_node = new Node(n_node,father_node);
						tree1.add_node(tmp_node);
						father_node++;
						check=0;
				}
				
				else if (line2[i]==':'){
					check=1;
					if (label != ""){
						n_leaf++;
						Leaf tmp_leaf = new Leaf(label,n_leaf,n_node);
						
					}
				}
				
				else if (line2[i]==','){
					check=0;
				}
				
				else if (line2[i]==')'){
					check=0;
					
				}
				
				else{
					if (check==0){
						label = label+line2[i];
					}
					else{
						b_length = b_length+line2[i];
					}
				}
			}
		}
		catch (IOException e) {
			e.printStackTrace();
		}
		
		return null;
	}
}
