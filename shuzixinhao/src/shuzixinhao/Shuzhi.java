package shuzixinhao;

import java.util.Scanner;

public class Shuzhi {

	public static void main(String[] args) {
		// TODO Auto-generated method stub
		DFT y=new DFT(16);
		
		long start=System.nanoTime();
		y.in();
		y.yunsaun();
		long end=System.nanoTime(); //获取结束时间
		System.out.println("DFT程序运行时间： "+(end-start)/1000000+"ms");
		FFT m=new FFT(16);
		
		long startTime=System.nanoTime();
		m.in();
		m.yunsuan();
		long endTime=System.nanoTime(); //获取结束时间
		System.out.println("FFT程序运行时间： "+(endTime-startTime)/1000000+"ms");
	}

}
class fushu{
	double shi,xu;
	public fushu(Double a,Double b) {
		// TODO Auto-generated constructor stub
		shi=a;
		xu=b;
	}
}
class caclu{
	fushu a[]=new fushu[2];
	Scanner input=new Scanner(System.in);
	void in(){
		for (int i = 0; i < a.length; i++) {
			a[i]=new fushu(input.nextDouble(), input.nextDouble());
		}
	}
	double[] jia(){
		double x[]=new double[2];
		x[0]=a[0].shi+a[1].shi;
		x[1]=a[0].xu+a[1].xu;
		return x;
		
	}
	double[] jian(){
		double x[]=new double[2];
		x[0]=a[0].shi-a[1].shi;
		x[1]=a[0].xu-a[1].xu;
		return x;
	}
	double[] cheng(){
		double x[]=new double[2];
		x[0]=a[0].shi*a[1].shi-a[0].xu*a[1].xu;
		x[1]=a[0].shi*a[1].xu+a[0].xu*a[1].shi;
		return x;
	}	
}
class FFT{
	fushu c[];
	fushu x[];
	FFT(int n){
		c=new fushu[n];
		x=new fushu[n];
	}
	void in(){
		for (int i = 0; i < c.length; i++) {
			c[i]=new fushu(i*1.0, 0.0);
		}
		for (int i = 0; i < x.length; i++) {
			x[i]=c[order(i)];
			
		}
		
	}
	int order(int n){//输入当前值与总点数
		int x=(int) (Math.log(c.length) / Math.log(2));
		int a[]=new int[x];
		int bin=0;
		for (int i = 0; i < x; i++) {
			a[i]=n%2;
			n=n/2;
		}
		for (int i = 0; i < x; i++) {
			if (a[i]>0) {
				bin=(int) (bin+Math.pow(2, x-1-i));	
			}
		}
		return bin;
	}
	fushu[] butter(fushu x,fushu y,fushu z){
		fushu m[]=new fushu[2];
		caclu nm=new caclu();
		nm.a[0]=new fushu(z.shi, z.xu);
		nm.a[1]=new fushu(y.shi, y.xu);
		fushu a=new fushu(nm.cheng()[0], nm.cheng()[1]);
		m[0]=new fushu(x.shi+a.shi, x.xu+a.xu);
		m[1]=new fushu(x.shi-a.shi, x.xu-a.xu);
		return m;
	}
	void yunsuan(){
		int M=(int) (Math.log(c.length) / Math.log(2));
		caclu mn=new caclu();
		for (int L = 1; L <= M; L++) {
			int B=(int)(Math.pow(2, L-1));
			for (int J = 0; J <= B-1; J++) {
				int P=(int)(Math.pow(2, M-L))*J;
				for (int k = J; k <= c.length-1; k+=(int)(Math.pow(2, L))) {
					fushu z=new fushu(Math.cos(2*Math.PI*P/c.length), -Math.sin(2*Math.PI*P/c.length));
					fushu y=new fushu(x[k+B].shi, x[k+B].xu);
					fushu n=new fushu(x[k].shi, x[k].xu);
					x[k]=butter(n, y, z)[0];
					x[k+B]=butter(n, y, z)[1];
				}
			}
		}
		for (int i = 0; i < x.length; i++) {
			System.out.println(x[i].shi+" "+x[i].xu);	
		}
	}
}
class DFT{
	double sum_shi=0;
	double sum_xu=0;
	fushu b[];
	public DFT(int n) {
	// TODO Auto-generated constructor stub
	b=new fushu[n];	
}
	void in(){
		for (int i = 0; i < b.length; i++) {
			b[i]=new fushu(i*1.0,0.0);
		}
	}
	void yunsaun(){
		caclu ab =new caclu();
		
		for (int i = 0; i < b.length; i++) {
			sum_shi=0;sum_xu=0;
			for (int j = 0; j < b.length; j++) {
				ab.a[1]=new fushu(Math.cos(2*Math.PI/b.length*i*j), -Math.sin(2*Math.PI/b.length*i*j));
				ab.a[0]=new fushu(b[j].shi, b[j].xu);
				sum_shi=sum_shi+ab.cheng()[0];
				sum_xu=sum_xu+ab.cheng()[1];
			}
			System.out.println(sum_shi+"  "+sum_xu);
		}
	}
}