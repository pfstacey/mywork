public class Node {
    protected String word;
    protected double count;

    public Node(String w, double c){
        this.word = w;
        this.count = c;
    }

    public String getWord(){
        return this.word;
    }

    public double getCount(){
        return this.count;
    }

    public void count(){
        ++this.count;
    }

    public void setCount(double c){
        this.count = c;
    }

    @Override
    public String toString(){
        return "[" + this.word + ", " + this.count + "]";
    }
}
