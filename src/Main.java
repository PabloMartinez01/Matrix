

public class Main {
    public static void main(String[] args) throws Exception {
       

        Matrix m = new Matrix(new double[][]{
            {1, -2, 2, 2},
            {0, 4, -2, 1},
            {1, -2, 4, 0},
            {1, -1, 2, 2}
        });

        Matrix m2 = new Matrix(new double[][]{
            {0.0, 1.0, 1.0, 5.0}, 
            {-2.0, 4.0, 2.0, 7.0},
            {1.0, 0.0, 2.0, 0.2},
            {2.0, 2.0, 3.0, 4.2}

        });

        System.out.println(m.add(m2));
        System.out.println(m.sub(m2));
        System.out.println(m.multiply(m2));
        System.out.println(m.det());
        System.out.println(m.rank());
        System.out.println(m.triangulate());
        System.out.println(m.inverse());
        System.out.println(m.multiply(m.inverse()));

        
    }

    
}
