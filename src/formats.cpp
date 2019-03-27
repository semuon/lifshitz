#include "formats.h"

void Formats::PrintVectorField(FILE *f, const double t, const VectorField &vec, const Lattice &lat)
{
 const uint count = vec.Count();
  uint dim= lat.Dim();
  const uint spdim = dim - pNumDiagDims;
  VECTOR<uint> L(dim);
  VECTOR<int> coor(dim);
  lat.LatticeSizes(L);
  uint dir=0; 

  for(uint i = 0; i < count ; i++)
  {
    lat.LinkCoordinates(coor, dir, i);
    if (dir==0)
    {
//       if (coor[0]==0 )
//       {
// 	SAFE_FPRINTF(f, "\n");
//       }
      SAFE_FPRINTF(f, "%2.4lf\t", t);
      for (uint j=0; j<spdim; j++)
      {
	SAFE_FPRINTF(f, "%d\t",coor[j]);
      }
    }
     if(dir==2)
     {
      SAFE_FPRINTF(f, "%2.15le\n", real(vec[i]));
     }
     else
     {
     SAFE_FPRINTF(f, "%2.15le\t", real(vec[i]));
     }
  }
  
//   SAFE_FPRINTF(f, "\n");

  fflush(f);
}

void Formats::PrintVector(FILE *f, const double t, const TBaseLinearVector<double> &vec)
{
  const uint count = vec.Count();

  SAFE_FPRINTF(f, "%2.4lf\t", t);

  for(uint i = 0; i < count - 1; i++)
  {
    SAFE_FPRINTF(f, "%2.15le\t", vec[i]);
  }
  SAFE_FPRINTF(f, "%2.15le\n", vec[count - 1]);

  fflush(f);
}

void Formats::PrintVector(FILE *f, const double t, const TBaseLinearVector<t_complex> &vec)
{
  const uint count = vec.Count();

  SAFE_FPRINTF(f, "%2.4lf\t", t);

  for(uint i = 0; i < count - 1; i++)
  {
    SAFE_FPRINTF(f, "%2.15le\t%2.15le\t", real(vec[i]), imag(vec[i]));
  }
  SAFE_FPRINTF(f, "%2.15le\t%2.15le\n", real(vec[count - 1]), imag(vec[count - 1]));

  fflush(f);
}

void Formats::PrintScalarField(FILE *f, const double t, const ScalarField  &val, const Lattice &lat)
{
  const uint count = val.Count();
  uint dim= lat.Dim();
  const uint spdim = dim - pNumDiagDims;
  VECTOR<uint> L(dim);
  VECTOR<int> coor(dim);
  lat.LatticeSizes(L);

  for(uint i = 0; i < count ; i++)
  {
    lat.SiteCoordinates(coor, i);
    if (coor[0]==0)
    {
      SAFE_FPRINTF(f, "\n");
    }
      SAFE_FPRINTF(f, "%2.4lf\t",t);
     for (uint j=0; j<spdim; j++)
     {
      SAFE_FPRINTF(f, "%d\t",coor[j]);
     }
     SAFE_FPRINTF(f, "  %2.15le \n", real(val[i]));
  }
  
  SAFE_FPRINTF(f, "\n");

  fflush(f);
}

void Formats::PrintScalar(FILE *f, const double t, const double val)
{
  SAFE_FPRINTF(f, "%2.4lf\t%2.15le\n", t, val);

  fflush(f);
}


void Formats::PrintScalar(FILE *f, const double t, const t_complex val)
{
  SAFE_FPRINTF(f, "%2.4lf\t%2.15le\t%2.15le\n", t, real(val), imag(val));

  fflush(f);
}


void Formats::ReadVectorField( const std::string name, const double t, VectorField &vec, const double rescale, const Lattice &lat)
{
 FILE *f=NULL; 
 const uint count = vec.Count();
 uint dim= lat.Dim();
 uint vindex=0; 
 const std::string f_attr = "r";
  const uint spdim = dim - pNumDiagDims;
 VECTOR<double> v_tmp(dim);
 VECTOR<int> coor_tmp(dim);
 double t_tmp=0.; 
 uint check=0;
 
 f = pDataDir.OpenFile(name, f_attr);
 
 std::cout<<"reading file:"<< name<< std::endl; 
 
 while(fgetc(f)!=EOF)
 {
   if (spdim == 0)
   {
    SAFE_FSCANF(f,"%lf %lf %lf %lf", &t_tmp, &v_tmp[0], &v_tmp[1], &v_tmp[2]);
   }
   else if (spdim == 1)
   {
    SAFE_FSCANF(f,"%lf  %d %lf %lf %lf", &t_tmp, &coor_tmp[0], &v_tmp[0], &v_tmp[1], &v_tmp[2]);
    std::cout <<" "<< t_tmp<<" " << coor_tmp[0]<<" " <<v_tmp[0]<<" " <<v_tmp[1]<<" " << v_tmp[2] << std::endl; 
     std::cout <<"check= "<<check << "count= "<<count << std::endl;
   }
   else if (spdim == 2)
   {
    SAFE_FSCANF(f,"%lf %d %d %lf %lf %lf", &t_tmp,&coor_tmp[0],&coor_tmp[1], &v_tmp[0], &v_tmp[1], &v_tmp[2]);
   }
   else if (spdim == 3)
   {
    SAFE_FSCANF(f,"%lf %d %d %d %lf %lf %lf", &t_tmp, &coor_tmp[0], &coor_tmp[1], &coor_tmp[2], &v_tmp[0], &v_tmp[1], &v_tmp[2]);
   }
   if (fabs(t_tmp -t) < 1.e-6)
   {
     for (uint dir=0; dir<dim; dir++)
     {
       vindex=lat.SiteIndex(coor_tmp);
       vec(vindex, dir) = v_tmp[dir];
       check++;
     }
   }
   if (t_tmp > t) break; 
 }

 std::cout <<"check= "<<check << "count= "<<count << std::endl;  
 ASSERTVERB(check==count,
	     "Formats::Reading vectors field from files\n \
	     the vectors field doesn't seem to have the right dimension\n");
 
 
 // rescaling
 vec*=rescale; 
 
 fclose(f); 
    
}

void Formats::PrintMatrix(FILE *f, const Matrix &m)
{
  const uint ncols = m.Ncols();
  const uint nrows = m.Nrows();

  for(uint i = 0; i < nrows; i++)
  {
    for(uint j = 0; j < ncols - 1; j++)
    {
      SAFE_FPRINTF(f, "%2.15le\t%2.15le\t", real(m(i, j)), imag(m(i, j)));
    }
    SAFE_FPRINTF(f, "%2.15le\t%2.15le\n", real(m(i, ncols - 1)), imag(m(i, ncols - 1)));
  }

  fflush(f);
}

