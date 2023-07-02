import java.math.BigDecimal;
import java.math.RoundingMode;
import java.util.Arrays;

public class Matrix {

    private double[][] data;
    private int rows;
    private int cols;

    public Matrix(double[][] data) {

        if (!isConsistent(data))
            throw new MatrixException("Data is empty or dimensions do not match");

        this.data = data;
        rows = data.length;
        cols = data[0].length;

    }


    public double get(int i, int j){

        if (!validRow(i) || !validCol(j))
            throw new MatrixException("Not valid i or j coord");

        return data[i][j];
    }


    public double[] getRow(int i){

        if (!validRow(i))
            throw new MatrixException("Not valid i coord");
        return data[i];
    }


    public double[] getCol(int j){

        if (!validCol(j))
            throw new MatrixException("Not valid j coord");

        double res[] = new double[data.length];
        for (int i = 0; i < data.length; i++){
            res[i] = data[i][j];
        }
        return res;
    }


    public Matrix add(Matrix matrix) {

        if (!equalDimensions(matrix))
            throw new MatrixException("Dimensions do not match");

        double res[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res[i][j] = data[i][j] + matrix.data[i][j];
            }
        }

        return new Matrix(res);
    }


    public Matrix sub(Matrix matrix) {

        if (!equalDimensions(matrix))
            throw new MatrixException("Dimensions do not match");

        return add(matrix.negate());
    }


    public Matrix multiply(double value) {

        double res[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                res[i][j] = data[i][j] * value;

        return new Matrix(res);

    }


    public Matrix multiply(Matrix matrix) {

        if (!isMultipliable(matrix))
            throw new MatrixException("Matrixes can not be multipied");

        double res[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {

                res[i][j] = 0;

                for (int k = 0; k < cols; k++)
                    res[i][j] += data[i][k] * matrix.data[k][j];
            }
        }

        return new Matrix(res).round(10);
    }


    public Matrix negate() {
        double res[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                res[i][j] = -data[i][j];

        return new Matrix(res);
    }


    public double det() {

        if (!isSquared())
            throw new MatrixException("Matrix is not squared");

        return det(this);
    }

    private double det(Matrix m) {

        if (m.data.length == 1) return m.data[0][0];
        if (m.data.length == 2) return m.data[0][0] * m.data[1][1] - m.data[0][1] * m.data[1][0];

        double total = 0;

        for (int j = 0; j < m.cols; j++) {

            double e = (j % 2 == 1) ? -1 : 1;
            double det = det(m.adj(0, j));

            total += e * m.data[0][j] * det;

        }

        return round(total, 10);
    }

    @Override
    public Matrix clone(){
        double m[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                m[i][j] = data[i][j];

        return new Matrix(m);
    }


    public Matrix triangulate() {

        double m[][] = clone().data;

        for (int j = 0; j < Math.min(m.length, m[0].length) ; j++) {

            if (m[j][j] == 0.0) 
                reposition(m, j);
                    
            for (int i = j + 1; i < m.length; i++) {

                if (m[i][j] != 0.0) {

                    double prop = m[j][j] / m[i][j];

                    for (int k = 0; k < m[i].length; k++) {
                        m[i][k] = m[i][k] - (m[j][k] / prop);
                    }
                }
            }
        }
        return new Matrix(m).round(10);
    }


    public Matrix transpose() {

        double res[][] = new double[cols][rows];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {
                res[j][i] = data[i][j];
            }
        }

        return new Matrix(res);
    }


    private int countNonZeroVectors() {
        int counter = 0;
        for (int i = 0; i < data.length; i++) {
            if (data[i][data[i].length - 1] != 0)
                counter++;
        }
        return counter;
    }
    

    public Matrix round(int precision) {

        double res[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++)
            for (int j = 0; j < cols; j++)
                res[i][j] = round(data[i][j], precision);

        return new Matrix(res);

    }

    public Matrix adj(int x, int y) {

        if (!validRow(x) || !validCol(y))
            throw new MatrixException("Not valid i or j coord");

        int ni = 0;
        int nj = 0;

        double res[][] = new double[rows - 1][cols - 1];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < rows; j++) {

                if (i != x && j != y) {
                    res[ni][nj] = data[i][j];
                    nj++;
                }
            }
            if (i != x)
                ni++;
            nj = 0;
        }

        return new Matrix(res);

    }


    public Matrix getAdjMatrix() {

        double res[][] = new double[rows][cols];

        for (int i = 0; i < rows; i++) {
            for (int j = 0; j < cols; j++) {

                int pow = ((i + j) % 2 == 0) ? 1 : -1;

                res[i][j] = pow * adj(i, j).det();

            }
        }

        return new Matrix(res);

    }


    public Matrix inverse() throws MatrixException {

        if (!isSquared())
            throw new MatrixException("Matrix is not squared");

        try {
            return getAdjMatrix().transpose().multiply(1 / det());
        } 
        catch (Exception e) {
            throw new MatrixException("Matrix is not invertible");
        }

    }

    public double rank() {
        return triangulate().countNonZeroVectors();
    }
    
    private boolean validRow(int i){
        return i >= 0 && i < data.length;
    }

    private boolean validCol(int j){
        return j >= 0 && j < data[0].length;
    }

    private boolean equalDimensions(Matrix matrix) {
        return rows == matrix.rows && cols == matrix.cols;
    }

    private boolean isMultipliable(Matrix matrix) {
        return cols == matrix.rows;
    }

    private boolean isSquared() {
        return rows == cols;
    }

    private static boolean isConsistent(double[][] data) {

        for (double[] row : data)
            if (row.length != data[0].length)
                return false;

        return data.length > 0;
    }

    private static void swap(double m[][], int i, int j){
        double aux[] = m[i];
        m[i] = m[j];
        m[j] = aux;
    }

     private static void reposition(double m[][], int i){

        for (int j = i + 1; j < m.length; j++){
            if (m[j][i] != 0.0){
                swap(m, i, j);
                break;
            }
        }
    }

    private static double round(double d, int precision) {
        return new BigDecimal(d).setScale(precision, RoundingMode.HALF_UP).doubleValue();
    }


    @Override
    public String toString() {

        StringBuilder stringBuilder = new StringBuilder();

        for (double[] r : data) {
            for (double c : r)
                stringBuilder.append(c + " ");
            stringBuilder.append("\n");
        }

        return new String(stringBuilder);
    }

    @Override
    public int hashCode() {
        final int prime = 31;
        int result = 1;
        result = prime * result + Arrays.deepHashCode(data);
        result = prime * result + rows;
        result = prime * result + cols;
        return result;
    }

    @Override
    public boolean equals(Object obj) {
        if (this == obj)
            return true;
        if (obj == null)
            return false;
        if (getClass() != obj.getClass())
            return false;
        Matrix other = (Matrix) obj;
        if (!Arrays.deepEquals(data, other.data))
            return false;
        if (rows != other.rows)
            return false;
        if (cols != other.cols)
            return false;
        return true;
    }

    

}
