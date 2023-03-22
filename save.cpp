#include "include/traj.h"
void save_files(int i_time, unsigned int n_space_div[3], float posL[3], float dd[3], double t,
                float np[2][n_space_divz][n_space_divy][n_space_divx], float currentj[2][3][n_space_divz][n_space_divy][n_space_divx],
                float V[n_space_divz][n_space_divy][n_space_divx],
                float E[3][n_space_divz][n_space_divy][n_space_divx], float B[3][n_space_divz][n_space_divy][n_space_divx],
                float KE[2][n_output_part], float posp[2][n_output_part][3])
{
#ifdef printDensity
  save_vti("Ne", i_time, n_space_div, posL, dd, n_cells, 1, t, (reinterpret_cast<const char *>(np[0])), "Float32", sizeof(float));
  // save_vti("Ne", i_time, n_space_div, posL, dd, n_cells, 1, t, np[0], "Float32", sizeof(float));
  save_vti_c("je", i_time, n_space_div, posL, dd, n_cells, 3, t, currentj[0], "Float32", sizeof(float));
#endif
#ifdef printV
  save_vti_c("V", i_time, n_space_div, posL, dd, n_cells, 1, t, V, "Float32", sizeof(float));
#endif
#ifdef printE
  save_vti_c("E", i_time, n_space_div, posL, dd, n_cells, 3, t, E, "Float32", sizeof(float));
#endif
#ifdef printB
  save_vti_c("B", i_time, n_space_div, posL, dd, n_cells, 3, t, B, "Float32", sizeof(float));
#endif
#ifdef printParticles
  save_vtp("e", i_time, n_output_part, 1, t, (reinterpret_cast<const char *>(&KE[0][0])), (reinterpret_cast<const char *>(&posp[0][0][0])));
  save_vtp("d", i_time, n_output_part, 1, t, (reinterpret_cast<const char *>(&KE[1][0])), (reinterpret_cast<const char *>(&posp[1][0][0])));
#endif
}
void save_vti(string filename, int i,
              unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
              const char *data, string typeofdata, int bytesperdata)
{
  std::ofstream os(outpath + filename + "_" + to_string(i) + ".vti", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n ";
  os << "<ImageData WholeExtent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\" ";
  os << "Origin=\"" + to_string(posl[0]) + " " + to_string(posl[1]) + " " + to_string(posl[2]) + "\"";
  os << " Spacing=\"" + to_string(dd[0]) + " " + to_string(dd[1]) + " " + to_string(dd[2]) + "\" ";
  os << "Direction=\"1 0 0 0 1 0 0 0 1\"> \n";
  os << "<FieldData>\n";
  os << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n";
  os << "</FieldData>\n";
  os << "<Piece Extent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\">\n";
  os << "<PointData Scalars=\"" + filename + "\">\n";
  os << "<DataArray type=\"" + typeofdata + "\" Name=\"" + filename + "\" NumberOfComponents=\"" + to_string(ncomponents) + "\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"16\" />\n";
  os << "  </PointData>\n";
  os << "<CellData>\n";
  os << "  </CellData>\n";
  os << "</Piece>\n";
  os << "</ImageData>\n";
  os << "<AppendedData encoding=\"raw\">_";
  uint64_t num1 = 8;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // single time double
  os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
  num1 = num * ncomponents * bytesperdata;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // data
  os.write(data, std::streamsize(num * ncomponents * bytesperdata));
  os << "</AppendedData>\n";
  os << "</VTKFile>";
  os.close();
}
/**
 * This corrects the order of dimensions for view in paraview, as opposed to save_vti which prints the raw data.
 */
void save_vti_c(string filename, int i,
                unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
                float data1[][n_space_divz][n_space_divy][n_space_divz], string typeofdata, int bytesperdata)
{
  if (ncomponents > 3)
  {
    cout << "Error: Cannot write file " << filename << " - too many components" << endl;
    return;
  }
  auto *data = new float[n_space_divz][n_space_divy][n_space_divx][3];
  for (int k = 0; k < n_space_div[2]; ++k)
  {
    for (int j = 0; j < n_space_div[1]; ++j)
    {
      for (int i = 0; i < n_space_div[0]; ++i)
      {
        for (int c = 0; c < ncomponents; ++c)
          data[k][j][i][c] = data1[c][k][j][i];
      }
    }
  }

  //  cout <<std::filesystem::temp_directory_path().string()+"/out/"+filename+"_"+to_string(i)+".vti";
  std::ofstream os(outpath + filename + "_" + to_string(i) + ".vti", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n ";
  os << "<ImageData WholeExtent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\" ";
  os << "Origin=\"" + to_string(posl[0]) + " " + to_string(posl[1]) + " " + to_string(posl[2]) + "\"";
  os << " Spacing=\"" + to_string(dd[0]) + " " + to_string(dd[1]) + " " + to_string(dd[2]) + "\" ";
  os << "Direction=\"1 0 0 0 1 0 0 0 1\"> \n";
  os << "<FieldData>\n";
  os << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n";
  os << "</FieldData>\n";
  os << "<Piece Extent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\">\n";
  os << "<PointData Scalars=\"" + filename + "\">\n";
  os << "<DataArray type=\"" + typeofdata + "\" Name=\"" + filename + "\" NumberOfComponents=\"" + to_string(ncomponents) + "\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"16\" />\n";
  os << "  </PointData>\n";
  os << "<CellData>\n";
  os << "  </CellData>\n";
  os << "</Piece>\n";
  os << "</ImageData>\n";
  os << "<AppendedData encoding=\"raw\">_";
  uint64_t num1 = 8;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // single time double
  os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
  num1 = num * ncomponents * bytesperdata;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // data
  os.write((reinterpret_cast<const char *>(data)), std::streamsize(num * ncomponents * bytesperdata));
  os << "</AppendedData>\n";
  os << "</VTKFile>";
  os.close();
  delete[] data;
}

void save_vti_c2(string filename, int i,
                 unsigned int n_space_div[3], float posl[3], float dd[3], uint64_t num, int ncomponents, double t,
                 float data1[3][n_space_divz2][n_space_divy2][n_space_divz2], string typeofdata, int bytesperdata)
{ // TODO: fix dimension sizes
  auto *data = new float[n_space_divz2][n_space_divy2][n_space_divx2][3];
  for (int k = 0; k < n_space_divz2; ++k)
  {
    for (int j = 0; j < n_space_divy2; ++j)
    {
      for (int i = 0; i < n_space_divx2; ++i)
      {
        for (int c = 0; c < 3; ++c)
          data[k][j][i][c] = data1[c][k][j][i];
      }
    }
  }

  //  cout <<std::filesystem::temp_directory_path().string()+"/out/"+filename+"_"+to_string(i)+".vti";
  std::ofstream os(outpath + filename + "_" + to_string(i) + ".vti", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"ImageData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\"> \n ";
  os << "<ImageData WholeExtent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\" ";
  os << "Origin=\"" + to_string(posl[0]) + " " + to_string(posl[1]) + " " + to_string(posl[2]) + "\"";
  os << " Spacing=\"" + to_string(dd[0]) + " " + to_string(dd[1]) + " " + to_string(dd[2]) + "\" ";
  os << "Direction=\"1 0 0 0 1 0 0 0 1\"> \n";
  os << "<FieldData>\n";
  os << "<DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n";
  os << "</FieldData>\n";
  os << "<Piece Extent=\"0 ";
  os << to_string(n_space_div[0] - 1) + " 0 " + to_string(n_space_div[1] - 1) + " 0 " + to_string(n_space_div[2] - 1) + "\">\n";
  os << "<PointData Scalars=\"" + filename + "\">\n";
  os << "<DataArray type=\"" + typeofdata + "\" Name=\"" + filename + "\" NumberOfComponents=\"" + to_string(ncomponents) + "\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"16\" />\n";
  os << "  </PointData>\n";
  os << "<CellData>\n";
  os << "  </CellData>\n";
  os << "</Piece>\n";
  os << "</ImageData>\n";
  os << "<AppendedData encoding=\"raw\">_";
  uint64_t num1 = 8;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // single time double
  os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
  num1 = num * ncomponents * bytesperdata;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  // data
  os.write((reinterpret_cast<const char *>(data)), std::streamsize(num * ncomponents * bytesperdata));
  os << "</AppendedData>\n";
  os << "</VTKFile>";
  os.close();
  delete[] data;
}
void save_vtp(string filename, int i, uint64_t num, int ncomponents, double t, const char *data, const char *points)
{
  std::ofstream os(outpath + filename + "_" + to_string(i) + ".vtp", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"PolyData\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  os << "<PolyData>\n <FieldData>\n  <DataArray type=\"Float64\" Name=\"TimeValue\" NumberOfTuples=\"1\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" offset=\"0\"/>\n</FieldData>\n";
  os << "<Piece NumberOfPoints=\"" + to_string(num);
  os << "\" NumberOfVerts=\"0\" NumberOfLines=\"0\" NumberOfStrips=\"0\" NumberOfPolys=\"0\" >\n";
  os << " <PointData>\n";
  os << "  <DataArray type=\"Float32\" Name=\" KE\" format=\"appended\" RangeMin=\"0\" RangeMax=\"0\" ";
  os << "offset=\"16\"/>\n";
  os << " </PointData>\n  <CellData>\n </CellData>\n";
  os << "<Points>";
  os << "    <DataArray type=\"Float32\" Name=\"Points\" NumberOfComponents=\"3\" format=\"appended\" RangeMin=\"0.0\" RangeMax=\"0.0\" offset=\"" + to_string(num * sizeof(float) * ncomponents + 24) + "\"/>\n";
  os << "  </Points>\n";
  os << "</Piece>\n";
  os << "</PolyData>\n";
  os << "<AppendedData encoding=\"raw\">\n_";
  uint64_t num1 = 8;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  os.write(reinterpret_cast<const char *>(&t), std::streamsize(sizeof(double)));
  num1 = num * ncomponents * sizeof(float);
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  os.write(data, std::streamsize(num * ncomponents * sizeof(float)));
  num1 = num * sizeof(float) * 3;
  os.write(reinterpret_cast<const char *>(&num1), std::streamsize(sizeof(num1)));
  os.write(points, std::streamsize(num * 3 * sizeof(float)));
  os << "\n</AppendedData>\n";
  os << "</VTKFile>";
  os.close();
}

void save_pvd(string filename, int ndatapoints)
{
  std::ofstream os(outpath + filename + ".pvd", std::ios::binary | std::ios::out);
  os << "<VTKFile type=\"Collection\" version=\"1.0\" byte_order=\"LittleEndian\" header_type=\"UInt64\">\n";
  os << "<Collection>\n";
  for (int i = 0; i < ndatapoints; i++)
  {
    os << "<DataSet timestep=\"" + to_string(i) + "\" part=\"0\" file=\"" + filename + "_" + to_string(i) + ".vti\"/>\n";
  }
  os << " </Collection>\n";
  os << "</VTKFile>\n";
  os.close();
}
