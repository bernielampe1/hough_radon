#include <exception>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string.h>
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

unsigned char getPixel(const int &numRows, const int &numCols,
                       const unsigned char *img, const int &r, const int &c) {
  if (r < numRows && r >= 0 && c < numCols && c >= 0) {
    return (img[r * numCols + c]);
  }

  return (0);
}

int integrateCircle(const int &x0, const int &y0, const int &r, const int &rows,
                    const int &cols, const unsigned char *img) {
  int f = 1 - r;
  int ddF_x = 1;
  int ddF_y = -2 * r;
  int x = 0;
  int y = r;
  int rp, cp;
  int integral = 0;

  rp = y0 + r;
  cp = x0;
  integral += getPixel(rows, cols, img, rp, cp);

  rp = y0 - r;
  cp = x0;
  integral += getPixel(rows, cols, img, rp, cp);

  rp = y0;
  cp = x0 + r;
  integral += getPixel(rows, cols, img, rp, cp);

  rp = y0;
  cp = x0 - r;
  integral += getPixel(rows, cols, img, rp, cp);

  while (x < y) {
    if (f >= 0) {
      y--;
      ddF_y += 2;
      f += ddF_y;
    }

    x++;
    ddF_x += 2;
    f += ddF_x;

    rp = y0 + y;
    cp = x0 + x;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 + y;
    cp = x0 - x;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 - y;
    cp = x0 + x;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 - y;
    cp = x0 - x;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 + x;
    cp = x0 + y;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 + x;
    cp = x0 - y;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 - x;
    cp = x0 + y;
    integral += getPixel(rows, cols, img, rp, cp);

    rp = y0 - x;
    cp = x0 - y;
    integral += getPixel(rows, cols, img, rp, cp);
  }

  return (integral);
}

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

class RadonSpace {
private:
  int **acc;
  int minRadius, maxRadius;
  int minCol, maxCol;
  int minRow, maxRow;
  int numRows, numCols, numRadii;

public:
  RadonSpace(const int minRad, const int maxRad, const int minC, const int maxC,
             const int minR, const int maxR)
      : minRadius(minRad), maxRadius(maxRad), minCol(minC), maxCol(maxC),
        minRow(minR), maxRow(maxR), numRows(maxRow - minRow + 1),
        numCols(maxCol - minCol + 1), numRadii(maxRadius - minRadius + 1) {
    acc = new int *[numRadii];

    for (int rad = 0; rad < numRadii; rad++) {
      acc[rad] = new int[numRows * numCols];
      memset((char *)(acc[rad]), '0', numRows * numCols * sizeof(int));
    }
  }

  void computeSpace(const int rows, const int cols, const unsigned char *img) {
    for (int r = 0; r < numRows; r++) {
      for (int c = 0; c < numCols; c++) {
        for (int rad = 0; rad < numRadii; rad++) {
          acc[rad][r * numCols + c] = integrateCircle(
              c + minCol, r + minRow, rad + minRadius, rows, cols, img);
        }
      }
    }
  }

  void writeSpace(const string &fname) {
    for (int rad = 0; rad < numRadii; rad++) {
      ofstream ofile;
      ostringstream oss;

      oss << fname << "_" << rad + minRadius << ".pgm";
      writePGM(oss.str(), numRows, numCols, acc[rad]);
    }
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

    RadonSpace acc(3, 50, 0, numCols - 1, 0, numRows - 1);

    acc.computeSpace(numRows, numCols, img);

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
