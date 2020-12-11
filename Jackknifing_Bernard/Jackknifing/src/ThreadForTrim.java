import java.util.ArrayList;
import java.util.HashMap;
import java.util.Random;


public class ThreadForTrim implements Runnable{

	private ArrayList<String> sequences;
	private int portion;
	private ArrayList<String> new_seq = new ArrayList<String> ();
	
	/**
	 * Thread constructor
	 * @param m_seq_posA position A in the list of seq.
	 * @param m_seq_posB position B in the list of seq.
	 * 
	 */
	public ThreadForTrim(ArrayList<String> m_sequences, int m_portion){
		this.sequences=m_sequences;
		this.portion=m_portion;
	}
	
	@Override
	public void run() {
		// TODO Auto-generated method stub
		for (int i=0;i<sequences.size();i++){
			int length = sequences.get(i).length();
			int new_length = (int)((length * portion) / 100);
			int del_size = length-new_length;
			
			StringBuilder tmp_seq = new StringBuilder(sequences.get(i));
			
			for (int j=0; j<del_size;j+=100){
				Random rand = new Random();
				int  random = rand.nextInt(length) + 1;
				if (random >= length-100){random = random - 100;}
				//tmp_seq = tmp_seq.deleteCharAt(random);
				int start = random;
				int end = random+100;
				tmp_seq.delete(start, end);
				length = tmp_seq.length();
			}
			this.new_seq.add(tmp_seq.toString());
		}
		
	}
	
	/**
	 * Return of list of k-mers library (k-mer library for each sequences).
	 * @return D2score list of k-mers library
	 * 
	 * @since 1.0
	 */
	public ArrayList<String> retNewSeq(){
		return new_seq;
	}

}
