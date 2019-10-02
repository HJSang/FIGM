
#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

//[[Rcpp::export]]
arma::mat mvrnormArma(int n, arma::vec mu, arma::mat sigma) {
  int ncols = sigma.n_cols;
  arma::mat Y = arma::randn(n, ncols);
  return arma::repmat(mu, 1, n).t() + Y * arma::chol(sigma);
}


//******************************************************
// Fractional Imputation using Gaussion Model
//********************************************************


//From Original data to extract obs and mis
void find_obs(mat & s,mat & obs, mat & mis)
{
  int n=s.n_rows;
  for(int i=0; i<n;i++)
  {
    uvec a=find((s.row(i)).t()==0);
    if(a.n_elem>0)
      mis=join_cols(mis,s.row(i));
    else
      obs=join_cols(obs,s.row(i));
  }
}



//Proposal disribution estimates
void proposal_est(mat & obs, vec & mu0, mat & S0)
{
  mu0=(mean(obs,0)).t();
  S0=cov(obs);
}


//Proposal Imputation
void Imputation(mat & s,vec & mu0, mat & S0, mat & com , int m)
{
  int n=s.n_rows;
  int i;
  vec z;
  uvec ind, ind1, ind2;
  mat temp;
  // add id variable
  mat s1=join_rows(linspace(1,n,n),s);
  for(i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    ind=find(z==0);
    if(ind.n_elem==0)
      temp=(join_cols(z, ones(2))).t();
    else
    {
      temp=repmat(z.t(),m,1);
      temp=join_rows(temp, zeros(m,2));
      temp.col(temp.n_cols-1)=1.0/m*ones(m);
      
      // imputation
      vec z1=z.subvec(1,z.n_elem-1);
      vec mu1=mu0(find(z1==0));
      vec mu2=mu0(find(z1!=0));
      int r1=mu1.n_elem;
      int r2=mu2.n_elem;
      ind1=find(z1==0);
      ind2=find(z1!=0);
      mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
      for(int i=0; i<r1; i++)
        for(int j=0; j<r1; j++)
          S1(i,j)=S0(ind1(i),ind1(j));
      
      for(int i=0; i<r2; i++)
        for(int j=0; j<r2; j++)
          S2(i,j)=S0(ind2(i),ind2(j));
      
      for(int i=0; i<r1; i++)
        for(int j=0; j<r2; j++)
          S3(i,j)=S0(ind1(i),ind2(j));
      vec b=z1(find(z1!=0));
      
      vec mu_star=mu1+S3*S2.i()*(b-mu2);
      mat S_star=S1-S3*S2.i()*S3.t();
      
      mat y_star=mvrnormArma(m,mu_star,S_star);
      for(int i=0; i<(int) ind.n_elem; i++)
        temp.col(ind(i))=y_star.col(i);
      
    }
    com=join_cols(com,temp);
    
  }
  
}

// parameter estimation

void est1(mat & com, vec & mu, mat & S, vec & prob, vec & q)
{
  int p=com.n_cols-3;
  vec w=com.col(com.n_cols-1);
  mat d=com.submat(0,1,com.n_rows-1,p);
  //cout<<"p\n"<<p<<"\n";
  for(int i=0;i<p; i++)
  {
    mu(i)=sum(d.col(i)%w)/sum(w);
    prob(i)=sum((d.col(i)<q(i))%w)/sum(w);
  }
  //cout<<"mu\n"<<mu<<"\n";
  
  for(int i=0;i<p; i++)
    for(int j=0; j<p; j++)
    {
      S(i,j)=sum((d.col(i)-mu(i))%(d.col(j)-mu(j))%w)/sum(w);
    }
    // cout<<"S\n"<<S<<"\n";
}



//weights update
void updatew(mat & com, vec & mu, mat & S, vec & mu0, mat & S0, int m)
{
  //int n=(com.col(0)).max();
  int p=com.n_cols;
  //cout<<"com\n"<<com.n_rows<<"\n";
  int i=0;
  double val0;
  double sign0;
  double val;
  double sign;
  
  log_det(val0, sign0, S0); 
  
  log_det(val, sign, S); 
  while(i <(int) com.n_rows)
  {
    if(com(i,p-2)==1)
    {
      i=i+1;
    }
    else
    {
      mat Y=com.submat(i,1,i+m-1,p-3);
      //cout<<"Y\n"<<Y<<"\n";
      vec w(m);
      for(int j=0; j<m; j++)
      {
        vec z=(Y.row(j)).t();
        w(j)=exp(0.5*as_scalar((z-mu0).t()*S0.i()*(z-mu0)-(z-mu).t()*S.i()*(z-mu))+
          0.5*val0-0.5*val);
      }
      w=w/sum(w);
      //cout<<"w\n"<<w<<"\n";
      com.submat(i,p-1,i+m-1,p-1)=w;
      i=i+m;
    }
    //cout<<"i\n"<<i<<"\n";
  }
}


//FI
//[[Rcpp::export]]
vec FIMV(mat s,vec q, int iter, int m)
{
  //int n=s.n_rows;
  int p=s.n_cols;
  int i;
  mat obs, mis;
  find_obs( s, obs,mis);
  vec mu0;
  mat S0;
  proposal_est( obs, mu0,S0);
  cout<<"mu0\t"<<mu0<<"\n";
  mat com;
  Imputation(s, mu0,S0,com ,m);
  //cout<<"com\n"<<com.head_rows(5)<<"\n";
  vec mu=mu0, prob=zeros(p);
  mat S(p,p);
  vec p1=mu0;
  
  for(i=0; i<iter; i++)
  {
    est1(com, mu, S, prob,q);
    cout<<"i\n"<<i<<"\n";
    cout<<"mu\n"<<mu<<"\n";
    updatew(com, mu, S, mu0, S0, m);
    if(norm(mu-p1)/norm(mu)<1e-6) break;
    p1=mu;
  }
  
  vec result=join_cols(mu, prob);
  result=join_cols(result, vectorise(S));
  
  return result;
}



//***********************************************
// Using GMM method
//************************************************



mat sub_mat(mat & com, vec  ID, int id)
{
  int p=com.n_cols,i;
  uvec ind=find(ID==id);
  int m=(int) ind.n_elem;
  mat temp(m,p);
  for(i=0; i<p; i++)
  {
    vec a=com.col(i);
    temp.col(i)=a(ind);
  }
  
  return temp;
  
}



//paramter estimates
void est(mat &com, mat & mu, mat & S, vec & prop, vec &mmu, vec & prob, vec & q)
{
  int p=S.n_cols;
  int G=prop.n_elem;
  mat e=com.submat(0,1,com.n_rows-1,p);
  for(int k=0; k<G; k++)
  {
    // cout<<"k\n"<<k<<"\n";
    mat tep=sub_mat(com,com.col(com.n_cols-2),k+1);
    vec w=tep.col(tep.n_cols-1);
    for(int j=0; j<p; j++)
      mu(k,j)=sum(tep.col(j+1)%w)/sum(w);
    prop(k)=sum(tep.col(tep.n_cols-1))/sum(com.col(com.n_cols-1));
  }
  
  for(int i=0; i<(int) e.n_rows; i++)
  {
    int g=com(i,com.n_cols-2);
    e.row(i)=e.row(i)-mu.row(g-1);
  }
  
  vec w1=com.col(com.n_cols-1);
  for(int i=0; i<p; i++)
    for(int j=0; j<p; j++)
    {
      S(i,j)=sum(e.col(i)%e.col(j)%w1)/sum(w1);
    }
    
    for(int i=0; i<(int) mmu.n_elem; i++)
    {
      mmu(i)=sum(com.col(i+1)%w1)/sum(w1);
      prob(i)=sum((com.col(i+1)<q(i))%w1)/sum(w1);
    }
    
    
}






vec BIC_cpp(mat & s, mat & mu, mat & S, vec & prop, int G )
{
  int n=s.n_rows;
  int p=s.n_cols;
  int i,k;
  double l=0.0;
  //cout<<"PG\t"<<PG<<"\n";
  for(i=0; i<n; i++)
  {
    vec y=(s.row(i)).t();
    //cout<<i<<"\t"<<y<<"\n";
    double f=0.0;
    uvec index=find(y==0);
    if(index.n_elem==0)
    {
      for(k=0; k<G; k++)
      {
        vec me=(mu.row(k)).t();
        mat va=S;
        f=f+prop(k)*pow(det(2*datum::pi*va),-0.5)*exp(-
          0.5*as_scalar((y-me).t()*va.i()*(y-me)));
      }
      
    }
    else
    {
      for(k=0; k<G; k++)
      {
        vec me=(mu.row(k)).t();
        vec me1=me(find(y!=0));
        mat va=S;
        uvec ind=find(y!=0);
        //cout<<"ind\t"<<ind<<"\n";
        mat va1(ind.n_elem, ind.n_elem);
        for(int u=0; u<(int)ind.n_elem; u++)
          for(int v=0; v<(int)ind.n_elem; v++)
            va1(u,v)=va(ind(u),ind(v));
        //cout<<"va1\t"<<va1<<"\n";
        //cout<<"va\t"<<va<<"\n";
        //cout<<"det\t"<<det(2*datum::pi*va1)<<"\n";
        f=f+prop(k)*pow(det(2*datum::pi*va1),-0.5)*exp(-
          0.5*as_scalar((y(ind)-me1).t()*va1.i()*(y(ind)-me1)));
      }
    }
    l=l+log(f);
  }
  
  vec cr(2);
  cr(0)=-2*l+log(n)*(G-1+G*p+p*(p+1)/2.0);
  cr(1)=-2*l+2*(G-1+G*p+p*(p+1)/2.0);
  return cr;
  
}


//Imputatin for model and Missing data
void update_ImputationGMM(mat & s, mat mu, mat S, int G, int m, mat & com)
{
  int n=s.n_rows;
  //int p=s.n_cols;
  mat s1=join_rows(linspace(1,n,n),s);
  //cout<<"s1\n"<<s1.head_rows(5)<<"\n";
  vec z;
  uvec ind;
  mat temp, com1;
  for(int i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    //cout<<"z\n"<<z<<"\n";
    ind=find(z==0);
    if(ind.n_elem==0)
    {
      temp=repmat(z.t(),G,1);
      temp=join_rows(temp,ones(G));
      temp=join_rows(temp,linspace(1,G,G));
      temp=join_rows(temp,1.0/G*ones(G));
      com1=join_cols(com1,temp);
    }
    else
    {
      for(int g=1; g<G+1; g++)
      {
        temp=repmat(z.t(),m,1);
        // cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,zeros(m));
        // cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,g*ones(m));
        //cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,1.0/(G*m)*ones(m));
        
        //cout<<"temp\n"<<temp<<"\n";
        // Imputation
        // imputation
        vec z1=z.subvec(1,z.n_elem-1);
        vec mu0=(mu.row(g-1)).t();
        vec mu1=mu0(find(z1==0));
        vec mu2=mu0(find(z1!=0));
        int r1=mu1.n_elem;
        int r2=mu2.n_elem;
        uvec ind1=find(z1==0);
        uvec ind2=find(z1!=0);
        mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
        for(int i=0; i<r1; i++)
          for(int j=0; j<r1; j++)
            S1(i,j)=S(ind1(i),ind1(j));
        
        for(int i=0; i<r2; i++)
          for(int j=0; j<r2; j++)
            S2(i,j)=S(ind2(i),ind2(j));
        
        for(int i=0; i<r1; i++)
          for(int j=0; j<r2; j++)
            S3(i,j)=S(ind1(i),ind2(j));
        vec b=z1(find(z1!=0));
        vec mu_star=mu1+S3*S2.i()*(b-mu2);
        mat S_star=S1-S3*S2.i()*S3.t();
        
        mat y_star=mvrnormArma(m,mu_star,S_star);
        //cout<<"y_star\n"<<y_star<<"\n";
        for(int i=0; i<(int) ind.n_elem; i++)
          temp.col(ind(i))=y_star.col(i);
        com1=join_cols(com1,temp);
        
      }
      
      
      
    }
    
  }
  com=com1;
}


void update_ImputationGMM2(mat & s, mat  mu, mat S, int G, int m, mat & com)
{
  int n=s.n_rows;
  //int p=s.n_cols;
  mat s1=join_rows(linspace(1,n,n),s);
  //cout<<"s1\n"<<s1.head_rows(5)<<"\n";
  vec z;
  uvec ind;
  mat temp, com1=com;
  int count=0;
  for(int i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    //cout<<"z\n"<<z<<"\n";
    ind=find(z==0);
    if(ind.n_elem==0)
    {
      count=count+G;
    }
    else
    {
      for(int g=1; g<G+1; g++)
      {
        temp=com1.submat(count,0,count+m-1,com1.n_cols-1);
        vec z1=z.subvec(1,z.n_elem-1);
        vec mu0=(mu.row(g-1)).t();
        vec mu1=mu0(find(z1==0));
        vec mu2=mu0(find(z1!=0));
        int r1=mu1.n_elem;
        int r2=mu2.n_elem;
        uvec ind1=find(z1==0);
        uvec ind2=find(z1!=0);
        mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
        for(int i=0; i<r1; i++)
          for(int j=0; j<r1; j++)
            S1(i,j)=S(ind1(i),ind1(j));
        
        for(int i=0; i<r2; i++)
          for(int j=0; j<r2; j++)
            S2(i,j)=S(ind2(i),ind2(j));
        
        for(int i=0; i<r1; i++)
          for(int j=0; j<r2; j++)
            S3(i,j)=S(ind1(i),ind2(j));
        vec b=z1(find(z1!=0));
        vec mu_star=mu1+S3*S2.i()*(b-mu2);
        mat S_star=S1-S3*S2.i()*S3.t();
        
        mat y_star=mvrnormArma(m,mu_star,S_star);
        //cout<<"y_star\n"<<y_star<<"\n";
        for(int i=0; i<(int) ind.n_elem; i++)
          temp.col(ind(i))=y_star.col(i);
        com1.submat(count,0,count+m-1,com1.n_cols-1)=temp;
        count=count+m;
        
      }
      
      
      
    }
    
  }
  com=com1;
}

//update weights
void updatew(mat & com, mat   mu, mat S, vec & prop, int G, int m)
{
  //int n=(com.col(0)).max();
  int p=com.n_cols;
  //cout<<"com\n"<<com.n_rows<<"\n";
  int i=0;
  double val;
  double sign;
  
  log_det(val, sign, S); 
  while(i <(int) com.n_rows)
  {
    if(com(i,p-3)==1)
    {
      mat Y=com.submat(i,1,i+G-1,p-4);
      
      vec w(G);
      
      for(int j=0; j<G; j++)
      {
        vec z=(Y.row(j)).t();
        vec mug= (mu.row(j)).t();
        w(j)=exp(-0.5*as_scalar((z-mug).t()*S.i()*(z-mug))-0.5*val)*prop(j);
      }
      w=w/sum(w);
      // cout<<"w\n"<<w<<"\n";
      com.submat(i,p-1,i+G-1,p-1)=w;
      
      i=i+G;
    }
    else
    {
      mat Y=com.submat(i,1,i+m*G-1,p-4);
      mat CV=cov(Y);
      uvec ind=find(CV.diag()>1e-8);
      vec z=(Y.row(0)).t();
      vec zobs=z(ind);
      mat S1=S(ind, ind);
      double val1;
      double sign1;
      
      log_det(val1, sign1, S1); 
      vec wg(G);
      for(int g=0; g<G; g++)
      {
        vec mug=(mu.row(g)).t();
        vec muobs=mug(ind);
        wg(g)=exp(-0.5*as_scalar(
          (zobs-muobs).t()*S1.i()*(zobs-muobs))-0.5*val1)*prop(g);
        
      }
      wg=wg/sum(wg);
      vec wf;
      for(int g=0; g<G; g++)
      {
        vec w(m);
        for(int j=0; j<m; j++)
        {
          w(j)=1.0/m;
        }
        w=w/sum(w)*wg(g);
        wf=join_cols(wf,w);
        
      }
      //cout<<"w\n"<<wf.head(5)<<"\n";
      com.submat(i,p-1,i+m*G-1,p-1)=wf;
      i=i+m*G;
    }
  }
}



//[[Rcpp::export]]
List FIGMM(mat & s,vec q, int iter, int m, int G, mat &mu0, mat & S0, vec & prop)
{
  //int n=s.n_rows;
  int p=s.n_cols;
  int i;
  mat obs, mis;
  find_obs( s, obs,mis);
  mat com;
  update_ImputationGMM(s,mu0,S0,G,m, com);
  mat mu=mu0;
  vec prob=zeros(p),mmu=(mean(mu0,0)).t();
  mat S=S0;
  vec p1=prop;
  wall_clock timer1, timer2, timer3;
  double t1=0,t2=0,t3=0;
  for(i=0; i<iter; i++)
  {
    
    timer1.tic();
    updatew(com, mu,S,prop,G, m);
    t1 = t1+timer1.toc();
    timer2.tic();
    est(com, mu,  S,  prop, mmu,  prob, q);
    t2 = t2+timer2.toc();
    //cout<<"i\n"<<i<<"\n";
    timer3.tic();
    update_ImputationGMM2(s,mu,S,G,m, com);
    t3 = t3+timer3.toc();
    if(norm(prop-p1)/norm(prop)<1e-6 &&i>100) break;
    p1=prop;
  }
  
  cout<<"time\n"<<t1<<"\n"<<t2<<"\n"<<t3<<"\n";
  vec bic=BIC_cpp( s,  mu, S,  prop, G );
  
  return List::create(Named("mu") = mu, Named("S")=S, Named("prop")=prop, Named("prob")=prob,
                            Named("overallmu")=mmu,Named("com")=com,Named("BIC")=bic);
  
}












// Imputation using GMM in the context of survey data

//Imputatin for model and Missing data
void update_ImputationGMMS(mat & s, mat & mu, mat & S, int G, int m, mat & com)
{
  int n=s.n_rows;
  //int p=s.n_cols;
  mat s1=join_rows(linspace(1,n,n),s);
  //cout<<"s1\n"<<s1.head_rows(5)<<"\n";
  vec z;
  uvec ind;
  mat temp, com1;
  for(int i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    //cout<<"z\n"<<z<<"\n";
    ind=find(z==0);
    if(ind.n_elem==0)
    {
      temp=repmat(z.t(),G,1);
      temp=join_rows(temp,ones(G));
      temp=join_rows(temp,linspace(1,G,G));
      temp=join_rows(temp,1.0/G*ones(G));
      com1=join_cols(com1,temp);
    }
    else
    {
      for(int g=1; g<G+1; g++)
      {
        temp=repmat(z.t(),m,1);
        // cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,zeros(m));
        // cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,g*ones(m));
        //cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,1.0/(G*m)*ones(m));
        
        //cout<<"temp\n"<<temp<<"\n";
        // Imputation
        // imputation
        vec z1=z.subvec(1,z.n_elem-2);
        vec mu0=(mu.row(g-1)).t();
        vec mu1=mu0(find(z1==0));
        vec mu2=mu0(find(z1!=0));
        int r1=mu1.n_elem;
        int r2=mu2.n_elem;
        uvec ind1=find(z1==0);
        uvec ind2=find(z1!=0);
        mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
        for(int i=0; i<r1; i++)
          for(int j=0; j<r1; j++)
            S1(i,j)=S(ind1(i),ind1(j));
        
        for(int i=0; i<r2; i++)
          for(int j=0; j<r2; j++)
            S2(i,j)=S(ind2(i),ind2(j));
        
        for(int i=0; i<r1; i++)
          for(int j=0; j<r2; j++)
            S3(i,j)=S(ind1(i),ind2(j));
        vec b=z1(find(z1!=0));
        vec mu_star=mu1+S3*S2.i()*(b-mu2);
        mat S_star=S1-S3*S2.i()*S3.t();
        
        mat y_star=mvrnormArma(m,mu_star,S_star);
        //cout<<"y_star\n"<<y_star<<"\n";
        for(int i=0; i<(int) ind.n_elem; i++)
          temp.col(ind(i))=y_star.col(i);
        com1=join_cols(com1,temp);
        
      }
      
      
      
    }
    
  }
  com=com1;
  com.col(com.n_cols-1)=com.col(com.n_cols-1)%com.col(com.n_cols-4);
}

//update weights
void updatewS(mat & com, mat &  mu, mat & S, vec & prop, int G, int m)
{
  //int n=(com.col(0)).max();
  int p=com.n_cols;
  //cout<<"com\n"<<com.n_rows<<"\n";
  int i=0;
  double val;
  double sign;
  
  log_det(val, sign, S); 
  while(i <(int) com.n_rows)
  {
    if(com(i,p-3)==1)
    {
      mat Y=com.submat(i,1,i+G-1,p-5);
      
      vec w(G);
      
      for(int j=0; j<G; j++)
      {
        vec z=(Y.row(j)).t();
        vec mug= (mu.row(j)).t();
        w(j)=exp(-0.5*as_scalar((z-mug).t()*S.i()*(z-mug))-0.5*val)*prop(j);
      }
      w=w/sum(w);
      // cout<<"w\n"<<w<<"\n";
      com.submat(i,p-1,i+G-1,p-1)=w;
      
      i=i+G;
    }
    else
    {
      mat Y=com.submat(i,1,i+m*G-1,p-5);
      mat CV=cov(Y);
      uvec ind=find(CV.diag()>1e-8);
      vec z=(Y.row(0)).t();
      vec zobs=z(ind);
      mat S1=S(ind, ind);
      double val1;
      double sign1;
      
      log_det(val1, sign1, S1); 
      vec wg(G);
      for(int g=0; g<G; g++)
      {
        vec mug=(mu.row(g)).t();
        vec muobs=mug(ind);
        wg(g)=exp(-0.5*as_scalar(
          (zobs-muobs).t()*S1.i()*(zobs-muobs))-0.5*val1)*prop(g);
        
      }
      wg=wg/sum(wg);
      vec wf;
      for(int g=0; g<G; g++)
      {
        vec w(m);
        for(int j=0; j<m; j++)
        {
          w(j)=1.0/m;
        }
        w=w/sum(w)*wg(g);
        wf=join_cols(wf,w);
        
      }
      //cout<<"w\n"<<wf.head(5)<<"\n";
      com.submat(i,p-1,i+m*G-1,p-1)=wf;
      i=i+m*G;
    }
  }
  com.col(com.n_cols-1)=com.col(com.n_cols-1)%com.col(com.n_cols-4);
}



void update_ImputationGMM2S(mat & s, mat & mu, mat & S, int G, int m, mat & com)
{
  int n=s.n_rows;
  //int p=s.n_cols;
  mat s1=join_rows(linspace(1,n,n),s);
  //cout<<"s1\n"<<s1.head_rows(5)<<"\n";
  vec z;
  uvec ind;
  mat temp, com1=com;
  int count=0;
  for(int i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    //cout<<"z\n"<<z<<"\n";
    ind=find(z==0);
    if(ind.n_elem==0)
    {
      count=count+G;
    }
    else
    {
      for(int g=1; g<G+1; g++)
      {
        temp=com1.submat(count,0,count+m-1,com1.n_cols-1);
        vec z1=z.subvec(1,z.n_elem-2);
        vec mu0=(mu.row(g-1)).t();
        vec mu1=mu0(find(z1==0));
        vec mu2=mu0(find(z1!=0));
        int r1=mu1.n_elem;
        int r2=mu2.n_elem;
        uvec ind1=find(z1==0);
        uvec ind2=find(z1!=0);
        mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
        for(int i=0; i<r1; i++)
          for(int j=0; j<r1; j++)
            S1(i,j)=S(ind1(i),ind1(j));
        
        for(int i=0; i<r2; i++)
          for(int j=0; j<r2; j++)
            S2(i,j)=S(ind2(i),ind2(j));
        
        for(int i=0; i<r1; i++)
          for(int j=0; j<r2; j++)
            S3(i,j)=S(ind1(i),ind2(j));
        vec b=z1(find(z1!=0));
        vec mu_star=mu1+S3*S2.i()*(b-mu2);
        mat S_star=S1-S3*S2.i()*S3.t();
        
        mat y_star=mvrnormArma(m,mu_star,S_star);
        //cout<<"y_star\n"<<y_star<<"\n";
        for(int i=0; i<(int) ind.n_elem; i++)
          temp.col(ind(i))=y_star.col(i);
        com1.submat(count,0,count+m-1,com1.n_cols-1)=temp;
        count=count+m;
        
      }
      
      
      
    }
    
  }
  com=com1;
  //com.col(com.n_cols-1)=com.col(com.n_cols-1)%com.col(com.n_cols-4);
}


vec BIC_cppS(mat & s, mat & mu, mat & S, vec & prop, int G )
{
  int n=s.n_rows;
  int p=s.n_cols-1;
  int i,k;
  double l=0.0;
  double N=0.0;
  //cout<<"PG\t"<<PG<<"\n";
  for(i=0; i<n; i++)
  {
    vec y=(s.row(i)).t();
    double w=y(y.n_elem-1);
    N=N+w;
    y=y.subvec(0, y.n_elem-2);
    //cout<<i<<"\t"<<y<<"\n";
    double f=0.0;
    uvec index=find(y==0);
    if(index.n_elem==0)
    {
      for(k=0; k<G; k++)
      {
        vec me=(mu.row(k)).t();
        mat va=S;
        f=f+prop(k)*pow(det(2*datum::pi*va),-0.5)*exp(-
          0.5*as_scalar((y-me).t()*va.i()*(y-me)));
      }
      
    }
    else
    {
      for(k=0; k<G; k++)
      {
        vec me=(mu.row(k)).t();
        vec me1=me(find(y!=0));
        mat va=S;
        uvec ind=find(y!=0);
        //cout<<"ind\t"<<ind<<"\n";
        mat va1(ind.n_elem, ind.n_elem);
        for(int u=0; u<(int)ind.n_elem; u++)
          for(int v=0; v<(int)ind.n_elem; v++)
            va1(u,v)=va(ind(u),ind(v));
        //cout<<"va1\t"<<va1<<"\n";
        //cout<<"va\t"<<va<<"\n";
        //cout<<"det\t"<<det(2*datum::pi*va1)<<"\n";
        f=f+prop(k)*pow(det(2*datum::pi*va1),-0.5)*exp(-
          0.5*as_scalar((y(ind)-me1).t()*va1.i()*(y(ind)-me1)));
      }
    }
    l=l+w*log(f);
  }
  
  vec cr(2);
  cr(0)=-2*l+log(N)*(G-1+G*p+p*(p+1)/2.0);
  cr(1)=-2*l+2*(G-1+G*p+p*(p+1)/2.0);
  return cr;
  
}

//[[Rcpp::export]]
List FIGMMS(mat & s,vec q, int iter, int m, int G, mat &mu0, mat & S0, vec & prop)
{
  //int n=s.n_rows;
  int p=s.n_cols-1;
  int i;
  mat obs, mis;
  find_obs( s, obs,mis);
  cout<<"obs\n"<<obs.head_rows(5)<<"\n";
  mat com;
  update_ImputationGMMS(s,mu0,S0,G,m, com);
  cout<<"com\n"<<com.head_rows(5)<<"\n";
  mat mu=mu0;
  vec prob=zeros(p),mmu=(mean(mu0,0)).t();
  mat S=S0;
  vec p1=prop;
  wall_clock timer1, timer2, timer3;
  double t1=0,t2=0,t3=0;
  for(i=0; i<iter; i++)
  {
    
    timer1.tic();
    updatewS(com, mu,S,prop,G, m);
    t1 = t1+timer1.toc();
    timer2.tic();
    est(com, mu,  S,  prop, mmu,  prob, q);
    t2 = t2+timer2.toc();
    cout<<"i\n"<<i<<"\n";
    timer3.tic();
    update_ImputationGMM2S(s,mu,S,G,m, com);
    t3 = t3+timer3.toc();
    if(norm(prop-p1)/norm(prop)<1e-6 &&i>100) break;
    p1=prop;
  }
  
  cout<<"time\n"<<t1<<"\n"<<t2<<"\n"<<t3<<"\n";
  vec bic=BIC_cppS( s,  mu, S,  prop, G );
  
  return List::create(Named("mu") = mu, Named("S")=S, Named("prop")=prop, Named("prob")=prob,
                            Named("overallmu")=mmu,Named("com")=com,Named("BIC")=bic);
  
}





// Imputation using GMM with categorical variables
//update weights
void updatew2(mat & com, mat &  mu, mat & S, mat & prop, int G, int m)
{
  //int n=(com.col(0)).max();
  int p=com.n_cols;
  //cout<<"com\n"<<com.n_rows<<"\n";
  int i=0;
  double val;
  double sign;
  
  log_det(val, sign, S); 
  while(i <(int) com.n_rows)
  {
    if(com(i,p-3)==1)
    {
      mat Y=com.submat(i,2,i+G-1,p-4);
      int indx=com(i,1);
      vec px=(prop.row(indx-1)).t();
      vec w(G);
      
      for(int j=0; j<G; j++)
      {
        vec z=(Y.row(j)).t();
        vec mug= (mu.row(j)).t();
        w(j)=exp(-0.5*as_scalar((z-mug).t()*S.i()*(z-mug))-0.5*val)*px(j);
      }
      w=w/sum(w);
      // cout<<"w\n"<<w<<"\n";
      com.submat(i,p-1,i+G-1,p-1)=w;
      
      i=i+G;
    }
    else
    {
      mat Y=com.submat(i,2,i+m*G-1,p-4);
      mat CV=cov(Y);
      uvec ind=find(CV.diag()>1e-8);
      vec z=(Y.row(0)).t();
      vec zobs=z(ind);
      
      int indx=com(i,1);
      vec px=(prop.row(indx-1)).t();
      
      mat S1=S(ind, ind);
      double val1;
      double sign1;
      
      log_det(val1, sign1, S1); 
      vec wg(G);
      for(int g=0; g<G; g++)
      {
        vec mug=(mu.row(g)).t();
        vec muobs=mug(ind);
        wg(g)=exp(-0.5*as_scalar(
          (zobs-muobs).t()*S1.i()*(zobs-muobs))-0.5*val1)*px(g);
        
      }
      wg=wg/sum(wg);
      vec wf;
      for(int g=0; g<G; g++)
      {
        vec w(m);
        for(int j=0; j<m; j++)
        {
          w(j)=1.0/m;
        }
        w=w/sum(w)*wg(g);
        wf=join_cols(wf,w);
        
      }
      //cout<<"w\n"<<wf.head(5)<<"\n";
      com.submat(i,p-1,i+m*G-1,p-1)=wf;
      i=i+m*G;
    }
  }
}


//paramter estimates
void est2(mat &com, mat & mu, mat & S, mat & prop, vec &mmu, vec & prob, vec & q)
{
  int p=S.n_cols;
  int G=prop.n_cols;
  mat e=com.submat(0,2,com.n_rows-1,p+1);
  vec ow=com.col(com.n_cols-1);
  for(int k=0; k<G; k++)
  {
    // cout<<"k\n"<<k<<"\n";
    mat tep=sub_mat(com,com.col(com.n_cols-2),k+1);
    vec w=tep.col(tep.n_cols-1);
    for(int j=0; j<p; j++)
      mu(k,j)=sum(tep.col(j+2)%w)/sum(w);
    for(int l=0; l<(int)prop.n_rows; l++)
      prop(l,k)=sum(w(find(tep.col(1)==l+1)))/sum(ow(find(com.col(1)==l+1)));
  }

  
  
  for(int i=0; i<(int) e.n_rows; i++)
  {
    int g=com(i,com.n_cols-2);
    e.row(i)=e.row(i)-mu.row(g-1);
  }
  
  vec w1=com.col(com.n_cols-1);
  for(int i=0; i<p; i++)
    for(int j=0; j<p; j++)
    {
      S(i,j)=sum(e.col(i)%e.col(j)%w1)/sum(w1);
    }
    
    for(int i=0; i<(int) mmu.n_elem; i++)
    {
      mmu(i)=sum(com.col(i+2)%w1)/sum(w1);
      prob(i)=sum((com.col(i+2)<q(i))%w1)/sum(w1);
    }
    
    
}

vec BIC_cpp2(mat & s, mat & mu, mat & S, mat & prop, int G )
{
  int n=s.n_rows;
  int p=s.n_cols;
  int i,k;
  double l=0.0;
  mat s1=s.cols(1,p-1);
  //cout<<"PG\t"<<PG<<"\n";
  for(i=0; i<n; i++)
  {
    vec y=(s1.row(i)).t();
    vec px=(prop.row(s(i,0)-1)).t();
    //cout<<i<<"\t"<<y<<"\n";
    double f=0.0;
    uvec index=find(y==0);
    if(index.n_elem==0)
    {
      for(k=0; k<G; k++)
      {
        vec me=(mu.row(k)).t();
        mat va=S;
        f=f+px(k)*pow(det(2*datum::pi*va),-0.5)*exp(-
          0.5*as_scalar((y-me).t()*va.i()*(y-me)));
      }
      
    }
    else
    {
      for(k=0; k<G; k++)
      {
        vec me=(mu.row(k)).t();
        vec me1=me(find(y!=0));
        mat va=S;
        uvec ind=find(y!=0);
        //cout<<"ind\t"<<ind<<"\n";
        mat va1(ind.n_elem, ind.n_elem);
        for(int u=0; u<(int)ind.n_elem; u++)
          for(int v=0; v<(int)ind.n_elem; v++)
            va1(u,v)=va(ind(u),ind(v));
        //cout<<"va1\t"<<va1<<"\n";
        //cout<<"va\t"<<va<<"\n";
        //cout<<"det\t"<<det(2*datum::pi*va1)<<"\n";
        f=f+px(k)*pow(det(2*datum::pi*va1),-0.5)*exp(-
          0.5*as_scalar((y(ind)-me1).t()*va1.i()*(y(ind)-me1)));
      }
    }
    l=l+log(f);
  }
  
  int nx=prop.n_rows;
  vec cr(2);
  cr(0)=-2*l+log(n)*((G-1)*nx+G*p+p*(p+1)/2.0);
  cr(1)=-2*l+2*((G-1)*nx+G*p+p*(p+1)/2.0);
  return cr;
  
}

//Imputatin for model and Missing data
void update_ImputationGMM3(mat & s, mat & mu, mat & S, int G, int m, mat & com)
{
  int n=s.n_rows;
  //int p=s.n_cols;
  mat s1=join_rows(linspace(1,n,n),s);
  //cout<<"s1\n"<<s1.head_rows(5)<<"\n";
  vec z;
  uvec ind;
  mat temp, com1;
  for(int i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    //cout<<"z\n"<<z<<"\n";
    ind=find(z==0);
    if(ind.n_elem==0)
    {
      temp=repmat(z.t(),G,1);
      temp=join_rows(temp,ones(G));
      temp=join_rows(temp,linspace(1,G,G));
      temp=join_rows(temp,1.0/G*ones(G));
      com1=join_cols(com1,temp);
    }
    else
    {
      for(int g=1; g<G+1; g++)
      {
        temp=repmat(z.t(),m,1);
        // cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,zeros(m));
        // cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,g*ones(m));
        //cout<<"temp\n"<<temp<<"\n";
        temp=join_rows(temp,1.0/(G*m)*ones(m));
        
        vec z1=z.subvec(2,z.n_elem-1);
        vec mu0=(mu.row(g-1)).t();
        vec mu1=mu0(find(z1==0));
        vec mu2=mu0(find(z1!=0));
        int r1=mu1.n_elem;
        int r2=mu2.n_elem;
        uvec ind1=find(z1==0);
        uvec ind2=find(z1!=0);
        mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
        for(int i=0; i<r1; i++)
          for(int j=0; j<r1; j++)
            S1(i,j)=S(ind1(i),ind1(j));
        
        for(int i=0; i<r2; i++)
          for(int j=0; j<r2; j++)
            S2(i,j)=S(ind2(i),ind2(j));
        
        for(int i=0; i<r1; i++)
          for(int j=0; j<r2; j++)
            S3(i,j)=S(ind1(i),ind2(j));
        vec b=z1(find(z1!=0));
        vec mu_star=mu1+S3*S2.i()*(b-mu2);
        mat S_star=S1-S3*S2.i()*S3.t();
        
        mat y_star=mvrnormArma(m,mu_star,S_star);
        //cout<<"y_star\n"<<y_star<<"\n";
        for(int i=0; i<(int) ind.n_elem; i++)
          temp.col(ind(i))=y_star.col(i);
        com1=join_cols(com1,temp);
        
      }
      
      
      
    }
    
  }
  com=com1;
}


void update_ImputationGMM4(mat & s, mat & mu, mat & S, int G, int m, mat & com)
{
  int n=s.n_rows;
  //int p=s.n_cols;
  mat s1=join_rows(linspace(1,n,n),s);
  //cout<<"s1\n"<<s1.head_rows(5)<<"\n";
  vec z;
  uvec ind;
  mat temp, com1=com;
  int count=0;
  for(int i=0; i<n; i++)
  {
    z=(s1.row(i)).t();
    //cout<<"z\n"<<z<<"\n";
    ind=find(z==0);
    if(ind.n_elem==0)
    {
      count=count+G;
    }
    else
    {
      for(int g=1; g<G+1; g++)
      {
        temp=com1.submat(count,0,count+m-1,com1.n_cols-1);
        vec z1=z.subvec(2,z.n_elem-1);
        vec mu0=(mu.row(g-1)).t();
        vec mu1=mu0(find(z1==0));
        vec mu2=mu0(find(z1!=0));
        int r1=mu1.n_elem;
        int r2=mu2.n_elem;
        uvec ind1=find(z1==0);
        uvec ind2=find(z1!=0);
        mat S1(r1,r1), S2(r2,r2), S3(r1,r2); 
        for(int i=0; i<r1; i++)
          for(int j=0; j<r1; j++)
            S1(i,j)=S(ind1(i),ind1(j));
        
        for(int i=0; i<r2; i++)
          for(int j=0; j<r2; j++)
            S2(i,j)=S(ind2(i),ind2(j));
        
        for(int i=0; i<r1; i++)
          for(int j=0; j<r2; j++)
            S3(i,j)=S(ind1(i),ind2(j));
        vec b=z1(find(z1!=0));
        vec mu_star=mu1+S3*S2.i()*(b-mu2);
        mat S_star=S1-S3*S2.i()*S3.t();
        
        mat y_star=mvrnormArma(m,mu_star,S_star);
        //cout<<"y_star\n"<<y_star<<"\n";
        for(int i=0; i<(int) ind.n_elem; i++)
          temp.col(ind(i))=y_star.col(i);
        com1.submat(count,0,count+m-1,com1.n_cols-1)=temp;
        count=count+m;
        
      }
      
      
      
    }
    
  }
  com=com1;
}

//[[Rcpp::export]]
List FIGMM2(mat & s,vec q, int iter, int m, int G, mat &mu0, mat & S0, mat & prop)
{
  //int n=s.n_rows;
  int i;
  mat obs, mis;
  find_obs( s, obs,mis);
  mat com;
  update_ImputationGMM3(s,mu0,S0,G,m, com);
  mat mu=mu0;
  vec prob=zeros(q.n_elem),mmu=(mean(mu0,0)).t();
  mat S=S0;
  mat p1=prop;
  wall_clock timer1, timer2, timer3;
  double t1=0,t2=0,t3=0;
  for(i=0; i<iter; i++)
  {
    
    timer1.tic();
    updatew2(com, mu,S,prop,G, m);
    t1 = t1+timer1.toc();
    timer2.tic();
    est2(com, mu,  S,  prop, mmu,  prob, q);
    t2 = t2+timer2.toc();
    timer3.tic();
    update_ImputationGMM4(s,mu,S,G,m, com);
    t3 = t3+timer3.toc();
    if(norm(prop-p1)/norm(prop)<1e-6 &&i>100) break;
    p1=prop;
    //cout<<"i\n"<<i<<"\n";
    //cout<<"prop\n"<<p1<<"\n";
  }
  
  cout<<"time\n"<<t1<<"\n"<<t2<<"\n"<<t3<<"\n";
  vec bic=BIC_cpp2( s,  mu, S,  prop, G );
  
  return List::create(Named("mu") = mu, Named("S")=S, Named("prop")=prop, Named("prob")=prob,
                            Named("overallmu")=mmu,Named("com")=com,Named("BIC")=bic);
  
}





