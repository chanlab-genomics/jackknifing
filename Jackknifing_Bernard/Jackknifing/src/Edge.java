import java.util.ArrayList;


public class Edge {
	private int value;
	private int father;
	private ArrayList<Integer>sons;
	
	public Edge(int m_value, int m_father, ArrayList<Integer>m_sons){
		this.value = m_value;
		this.father = m_father;
		this.sons = m_sons;
	}//constructor
}
