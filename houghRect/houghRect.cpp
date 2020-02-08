#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

using namespace std;

class Exception : public exception {
private:
  string errStr;

public:
  Exception(const string &mes) : errStr(mes) {}

  ~Exception() throw() {}

  virtual const char *what() const throw() { return (this->errStr.c_str()); }
};

void readPGM(const string &fname, int &numRows, int &numCols,
             unsigned char **img) {
  ifstream ifile;
  int maxVal;
  char magic[2];

  ifile.open(fname.c_str());
  if (!ifile) {
    throw Exception("could not open input image for reading");
  }

  ifile.read(magic, 2 * sizeof(unsigned char));
  if (magic[0] != 'P' || magic[1] != '5') {
    throw Exception("input image not a pgm (P5)");
  }

  ifile >> numCols >> numRows >> maxVal;
  *img = new unsigned char[numCols * numRows];

  ifile.read((char *)(*img), numRows * numCols * sizeof(unsigned char));
}

void writePGM(const string &fname, const int &numRows, const int &numCols,
              const int *img) {
  int minVal, maxVal;
  unsigned char val;

  minVal = maxVal = img[0];
  for (int i = 0; i < numRows * numCols; i++) {
    if (img[i] < minVal)
      minVal = img[i];
    if (img[i] > maxVal)
      maxVal = img[i];
  }

  double scale = 255.0 / (maxVal > minVal ? (maxVal - minVal) : 1);

  ofstream ofile;
  ofile.open(fname.c_str());
  if (!ofile) {
    throw Exception("could not open file for writing");
  }

  ofile << "P5\n" << numCols << " " << numRows << "\n255\n";

  for (int i = 0; i < numRows * numCols; i++) {
    val = (unsigned char)((img[i] - minVal) * scale + 0.5);
    ofile.write((char *)&val, sizeof(unsigned char));
  }
}

class HoughSpace {
private:
  int ***acc;
  int minCol, maxCol;
  int minRow, maxRow;
  int minWidth, maxWidth;
  int minHeight, maxHeight;
  int numWidths, numHeights;
  int numRows, numCols;

public:
  HoughSpace(const int minW, const int maxW, const int minH, const int maxH,
             const int minC, const int maxC, const int minR, const int maxR)
      : acc(0), minCol(minC), maxCol(maxC), minRow(minR), maxRow(maxR),
        minWidth(minW), maxWidth(maxW), minHeight(minH), maxHeight(maxH),
        numWidths(maxW - minW + 1), numHeights(maxH - minH + 1),
        numRows(maxRow - minRow + 1), numCols(maxCol - minCol + 1) {
    acc = new int **[numWidths];

    for (int w = 0; w < numWidths; w++) {
      acc[w] = new int *[numHeights];

      for (int h = 0; h < numHeights; h++) {
        acc[w][h] = new int[numRows * numCols];
        memset((char *)(acc[w][h]), '0', numRows * numCols * sizeof(int));
      }
    }
  }

  ~HoughSpace() {
    for (int w = 0; w < numWidths; w++) {
      for (int h = 0; h < numHeights; h++) {
        delete[] acc[w][h];
      }

      delete[] acc[w];
    }

    delete[] acc;
  }

  void addVotes(const int inr, const int inc) {
    int h_2, w_2;
    int rp, cp;
    int rp_1, rp_2;
    int cp_1, cp_2;

    for (int w = 0; w < numWidths; w++) {
      w_2 = w >> 1;

      for (int h = 0; h < numHeights; h++) {
        h_2 = h >> 1;

        // for this height and width, vote for each rectangle center
        cp_1 = inc - w_2;
        cp_2 = inc + w_2;
        for (int r = 0; r <= h; r++) {
          rp = inr + r - h_2;

          if (rp >= 0 && rp < numRows && cp_1 >= 0 && cp_1 < numCols &&
              cp_2 >= 0 && cp_2 < numCols) {
            acc[w][h][rp * numCols + cp_1]++;
            acc[w][h][rp * numCols + cp_2]++;
          }
        }

        rp_1 = inr - h_2;
        rp_2 = inr + h_2;
        for (int c = 0; c <= w; c++) {
          cp = inc + c - w_2;

          if (cp >= 0 && cp < numCols && rp_1 >= 0 && rp_1 < numRows &&
              rp_2 >= 0 && rp_2 < numRows) {
            acc[w][h][rp_1 * numCols + cp]++;
            acc[w][h][rp_2 * numCols + cp]++;
          }
        }
      }
    }
  }

  void writeSpace(const string &fname) {
    // write projection of (w, h) , (w, r), (w, c), (h, r), (h, c), (r, c)
    ostringstream oss;
    int val;

    int *w_h_img = new int[numWidths * numHeights];
    int *w_r_img = new int[numWidths * numRows];
    int *w_c_img = new int[numWidths * numCols];
    int *h_r_img = new int[numHeights * numRows];
    int *h_c_img = new int[numHeights * numCols];
    int *r_c_img = new int[numRows * numCols];

    memset(w_h_img, '0', numWidths * numHeights * sizeof(int));
    memset(w_r_img, '0', numWidths * numRows * sizeof(int));
    memset(w_c_img, '0', numWidths * numCols * sizeof(int));
    memset(h_r_img, '0', numHeights * numRows * sizeof(int));
    memset(h_c_img, '0', numHeights * numCols * sizeof(int));
    memset(r_c_img, '0', numRows * numCols * sizeof(int));

    for (int w = 0; w < numWidths; w++) {
      for (int h = 0; h < numHeights; h++) {
        for (int r = 0; r < numRows; r++) {
          for (int c = 0; c < numCols; c++) {
            val = acc[w][h][r * numCols + c];

            w_h_img[w * numHeights + h] += val;
            w_r_img[w * numRows + r] += val;
            w_c_img[w * numCols + c] += val;
            h_r_img[h * numRows + r] += val;
            h_c_img[h * numCols + c] += val;
            r_c_img[r * numCols + c] += val;
          }
        }
      }
    }

    oss << fname << "_w_h.pgm";
    writePGM(oss.str(), numWidths, numHeights, w_h_img);

    oss.clear();
    oss.str("");
    oss << fname << "_w_r.pgm";
    writePGM(oss.str(), numWidths, numRows, w_r_img);

    oss.clear();
    oss.str("");
    oss << fname << "_w_c.pgm";
    writePGM(oss.str(), numWidths, numCols, w_c_img);

    oss.clear();
    oss.str("");
    oss << fname << "_h_r.pgm";
    writePGM(oss.str(), numHeights, numRows, h_r_img);

    oss.clear();
    oss.str("");
    oss << fname << "_h_c.pgm";
    writePGM(oss.str(), numHeights, numCols, h_c_img);

    oss.clear();
    oss.str("");
    oss << fname << "_r_c.pgm";
    writePGM(oss.str(), numRows, numCols, r_c_img);

    delete[] w_h_img;
    delete[] w_r_img;
    delete[] w_c_img;
    delete[] h_r_img;
    delete[] h_c_img;
    delete[] r_c_img;
  }
};

int main(int argc, char **argv) {
  unsigned char *img;
  int numRows, numCols;

  if (argc != 2) {
    cerr << "Usage: " << argv[0] << " <img.pgm>" << endl;
    return (1);
  }

  try {
    readPGM(argv[1], numRows, numCols, &img);

    HoughSpace acc(5, 100, 5, 100, 0, numCols - 1, 0, numRows - 1);

    for (int r = 0; r < numRows; r++) {
      for (int c = 0; c < numCols; c++) {
        if (img[r * numCols + c]) {
          acc.addVotes(r, c);
        }
      }
    }

    acc.writeSpace("accumulators");
  } catch (Exception &e) {
    cerr << e.what() << endl;
    return (1);
  } catch (bad_alloc) {
    cerr << "caught a bad allocation...out of memory" << endl;
    return (1);
  } catch (...) {
    cerr << "caught unhandled exception" << endl;
    return (1);
  }

  return (0);
}
