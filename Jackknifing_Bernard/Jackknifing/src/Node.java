import java.util.ArrayList;

public class Node {
	
	private int number;
	private int father;
	private ArrayList<Integer> sons;
	
	public Node(int m_number, int m_father){
		this.number = m_number;
		this.father = m_father;
	}//constructor
}
