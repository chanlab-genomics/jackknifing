import java.util.ArrayList;
import java.util.List;
import java.util.Random;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.TimeUnit;


public class JackKnifing {
	
	private ArrayList<ArrayList<String>> tmp_matrix = new ArrayList<ArrayList<String>>();
	private int del_number;
	
	/**
	 * List that contains sequences.
	 */
	static ArrayList<String> sequences = new ArrayList<String>();
	
	/**
	 * Number of thread used to complete the D2scoreMatrix
	 */
	static int nbr_thread=Runtime.getRuntime().availableProcessors();
	
	public JackKnifing(ArrayList<ArrayList<String>> m_matrix, Integer m_del){
		
		//tmp_matrix=(ArrayList<ArrayList<String>>) m_matrix.clone();
		for (int i=0;i<m_matrix.size();i++){
			tmp_matrix.add((ArrayList<String>) m_matrix.get(i).clone());
			System.out.println(m_matrix.get(i).size());
		}
		
		del_number = m_del;
	}
	
	
	public JackKnifing(ArrayList<String> m_sequences){
			
			for (int i=0;i<m_sequences.size();i++){
				sequences.add((String) m_sequences.get(i));
			}
			
	}
	
	
	
	public ArrayList<ArrayList<String>> matrixsampling(Integer n_row){
		
		//Random rand = new Random();
		int N =tmp_matrix.size();
		//int  random = rand.nextInt(N) + 1;
		//tmp_matrix.remove(random-1);
		tmp_matrix.remove(n_row-1);
		N = tmp_matrix.size();
		for (int i=0; i<N;i++){
			if (tmp_matrix.get(i).size() >= n_row){
				tmp_matrix.get(i).remove(n_row);
			}
		}
		
		return tmp_matrix;
		
	}
	
	public ArrayList<String> datasampling(Integer m_portion){
			
			
			int N =sequences.size();
			int portion = m_portion;
			ArrayList<String> new_seqs = new ArrayList<String> ();
			
			List<Runnable> runnables = new ArrayList<Runnable>();
			
			int nbr_seq = N/nbr_thread;
			
			int count = 0;
			for (int i=0; i<nbr_thread;i++){
				ArrayList<String> tmp_seq = new ArrayList<String>();
				for (int k=0;k<nbr_seq;k++){
					System.out.println("normal length: "+sequences.get(count).length());
					tmp_seq.add(sequences.get(count));
					count++;
				}
				runnables.add(new ThreadForTrim(tmp_seq, portion));
			}
			if ((N-count)!=0){
				ArrayList<String> tmp_seq = new ArrayList<String>();
				for (int i=count; i<N;i++){
					System.out.println("normal length: "+sequences.get(count).length());
					tmp_seq.add(sequences.get(count));
					count++;
				}
				runnables.add(new ThreadForTrim(tmp_seq, portion));
			}
			nbr_thread=runnables.size();
			ExecutorService execute = Executors.newFixedThreadPool(nbr_thread);
			executeThreads(execute, runnables);
			
			for (int j=0;j<nbr_thread;j++){
				ArrayList<String>tmp_new_seq = ((ThreadForTrim) runnables.get(j)).retNewSeq();
				for (int k=0;k<tmp_new_seq.size();k++){
					new_seqs.add(tmp_new_seq.get(k));
					System.out.println("new length: "+tmp_new_seq.get(k).length());
				}
			}
			
			System.out.println("new_data: "+new_seqs.size());
			
			return new_seqs;
			
	}
	
	/**
	 * Function used to execute a "pool" of thread (limit max fix to the number of cpu).
	 * @param service
	 * @param runnables List of threads.
	 * 
	 * @since 1.0
	 */
	public static void executeThreads(final ExecutorService service, List<Runnable> runnables){
		
		for(Runnable r : runnables){
			service.execute(r);
		}
		
		service.shutdown();
		try {
			service.awaitTermination(2, TimeUnit.DAYS);
		} catch (InterruptedException e) {
			e.printStackTrace();
		}
	}
}
