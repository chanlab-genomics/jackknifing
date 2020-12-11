import java.util.ArrayList;


public class Tree {
	private ArrayList<Node> nodes;
	private ArrayList<Leaf> leafs;
	private ArrayList<Edge> edges;
	
	public Tree(){
	}//constructor
	
	public void add_node(Node m_node){
		nodes.add(m_node);
	}
}
