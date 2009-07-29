


/* applyMorphAnd applies a morphological operation using the given mask and kernel. The suffix "And"
 * implies that the pixel at (i, j) is set if and only if every bit in the kernel that falls in the
 * support of the mask, matches every corresponding bit in the neigborhood of (i, j) in the image. */
template<class T>
void applyMorphAnd(T &input, binaryImage3d &morphkernel, binaryImage3d &mask, T &output) {
  int irows = input.rows();
  int icols = input.cols();
  int ilays = input.layers();

  output.setRows(irows);
  output.setCols(icols);
  output.setLayers(ilays);
  output.resize();

  int mrows = morphkernel.rows();
  int mcols = morphkernel.cols();
  int mlays = morphkernel.layers();
  int mrowsOver2 = mrows/2;
  int mcolsOver2 = mcols/2;
  int mlaysOver2 = mlays/2;
	
  for (int l = 0; l < ilays; l++) {
    int ls = l - mlaysOver2;
	for (int i = 0; i < irows; i++) {
      int is = i - mrowsOver2;
      for (int j = 0; j < icols; j++) {
        int js = j - mcolsOver2;

        bool val = 1;
        for (int ll = max(0, -ls); ll < min(mlays, ilays - ls); ll++) {
          for (int ii = max(0, -is); ii < min(mrows, irows - is); ii++) {
            for (int jj = max(0, -js); jj < min(mcols, icols - js); jj++) {
              if (mask.getVoxel(ii, jj, ll)) {
                val &= (input.getVoxel(is+ii, js+jj, ls+ll) == morphkernel.getVoxel(ii, jj, ll));
              }
            }
          }
        }
        output.setVoxel(i, j, l, val);	
      }
	}
  }
  return;
}

