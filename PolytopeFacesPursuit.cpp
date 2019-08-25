#include "PolytopeFacesPursuit.h"


void PolytopeFacesPursuit(VectorXd& hyb_observe, VectorXd& w, VectorXd& ScaleFactor, BOOLEAN IsHybrid)
{

ifstream ifs;
if (IsHybrid)
  ifs.open(KERNEL_MATRIX_H,ios::in | ios::binary );
else
  ifs.open(KERNEL_MATRIX_B,ios::in | ios::binary );
  
if (ifs.is_open()){
static double STemp[SEQ_NUMBER*SEQ_NUMBER];

ifs.read( (char *) &STemp, sizeof STemp);
ifs.close();
Map<MatrixXd> S(STemp,SEQ_NUMBER,SEQ_NUMBER);

Matrix2d u_sub;

// INITIALIZATION
int num_spe = hyb_observe.size(),maxIters = 10*num_spe, counter=0, Picked=0, reject_num=1, recruits=1, negweith=-1, first_spe, projectsetlocate, i,k;
double Threshold = 1e-1, norm_diag, maxScore, curscore;

hyb_observe.maxCoeff(&first_spe);
VectorXi Selected = VectorXi::Constant(num_spe,-1), Rejected = VectorXi::Constant(num_spe,-1);
VectorXd Ux=VectorXd::Zero(num_spe), hyb_vec=VectorXd::Zero(num_spe), x=VectorXd::Zero(num_spe), z=VectorXd::Zero(num_spe), polyface=VectorXd::Ones(num_spe), Spivot=VectorXd::Ones(num_spe), r_proj(num_spe), searchfaces(num_spe), polybound=VectorXd::Ones(num_spe), decomp=VectorXd::Zero(num_spe);
MatrixXd U = MatrixXd::Zero(num_spe,num_spe), projectset = MatrixXd::Zero(num_spe,num_spe);
projectset.col(0) = S.col(first_spe);
Selected(0)=first_spe;
hyb_vec(0)=hyb_observe(first_spe);
U(0,0) = sqrt(S(first_spe,first_spe));

//.............................................


while(counter<maxIters )//
{
  while(1)
  {
     x.head(recruits)=U.topLeftCorner(recruits,recruits).transpose().triangularView<Lower>().solve(hyb_vec.head(recruits));  
     U.topLeftCorner(recruits,recruits).triangularView<Upper>().solveInPlace(x.head(recruits));//x.head(recruits)=
     
     negweith=-1;
     for(i=0;i<recruits;i++)
       {if(x(i)<0)
         {negweith=i;
         break;}}
     if(negweith !=-1)
           { if(negweith==recruits-1)
               { Rejected(reject_num-1) = Selected(recruits-1);
                reject_num++;}
            else
               { reject_num=1;
            	U.block(0,negweith,recruits,recruits-negweith-1).noalias()
 = U.block(0,negweith+1,recruits,recruits-negweith-1);
            	projectset.block(0,negweith,num_spe,recruits-negweith-1).noalias()
 = projectset.block(0,negweith+1,num_spe,recruits-negweith-1);
            	hyb_vec.segment(negweith,recruits-negweith-1).noalias()
 = hyb_vec.segment(negweith+1,recruits-negweith-1);
            	Selected.segment(negweith,recruits-negweith+1).noalias()
 = Selected.segment(negweith+1,recruits-negweith+1);
            	//
               // cout << Selected.head(recruits)<< endl;
               // cout << endl;
                //
            	
            	for(k=negweith;k<recruits-1;k++)
                  { if(U(k+1,k)>0)
                      {norm_diag= U.block(k,k,2,1).norm();
                      //
                      //cout << norm_diag<< endl;
                      //
                       if(k<recruits-2)
                         {u_sub(0,0)=U(k,k); u_sub(1,1)=U(k,k); u_sub(0,1)=U(k+1,k); u_sub(1,0)=-U(k+1,k);
                         U.block(k,k+1,2,recruits-k-2) = (u_sub/norm_diag)*U.block(k,k+1,2,recruits-k-2);//-2
                         }
                       U(k,k)=norm_diag; U(k+1,k)=0;
                      }
                       else
                         continue;
                }
                
                }
            recruits--;}
        else
           { if(reject_num>1)
               {for(k=0;k<reject_num;k++)
                    Rejected(i)=-1;
                reject_num=1;
              }
            break;
           }
 
  }
  
 
  
  r_proj.noalias()
 = hyb_observe - projectset.leftCols(recruits)*x.head(recruits);
  //r_proj = hyb_observe - projectset.topLeftCorner(num_spe,recruits)*x.head(recruits);
  
  
  if(r_proj.maxCoeff()>Threshold)
     {maxScore=-10E5;
      searchfaces.noalias()
 = polybound-projectset.leftCols(recruits)*z.head(recruits);
     for(i=0;i<recruits;i++)
         {if(Selected(i)>=0)
            r_proj(Selected(i)) = -1;
         }
         
     for(i=0;i<reject_num;i++)
         {if(Rejected(i)>=0)
            r_proj(Rejected(i)) = -1;
         }
     for(i=0;i<num_spe;i++)   
         {if(r_proj(i)>0)
            {curscore = r_proj(i)/searchfaces(i);
             if(curscore>maxScore)
               {maxScore=curscore;
                Picked=i;

               }
             }
         }
     }
    else
        break;
        
    
    //
   // cout << Picked<< endl;
                      //

    projectsetlocate = 0;
    for(i=0;i<num_spe;i++)
        {if(Selected(i)>=0)
          {Spivot(projectsetlocate) = S(Selected(i),Picked);
           projectsetlocate++;
          }
        }
    
    decomp.head(recruits) = U.topLeftCorner(recruits,recruits).transpose().triangularView<Lower>().solve(Spivot.head(recruits));
    U.block(0,recruits,recruits,1)= decomp.head(recruits);
    U(recruits,recruits) = sqrt(S(Picked,Picked)- decomp.head(recruits).squaredNorm());
    hyb_vec(recruits) = hyb_observe(Picked);
    projectset.col(recruits)= S.col(Picked);
    Selected(recruits) =  Picked;
    recruits++;
    counter++;

}

for(i=0;i<recruits;i++)  
   w(Selected(i))= x(i)/ScaleFactor(Selected(i));
 w/=w.sum();
  
 }
  else
  cout << "Database file missing!"<< endl;

}

