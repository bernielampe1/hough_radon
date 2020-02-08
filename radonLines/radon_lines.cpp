#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>

using namespace std;

class Exception : public exception {
private:
  const char *message;

public:
  Exception(const char *m) : message(m) {}

  virtual const char *what() const throw() { return (this->message); }
};

class PGMImage {
private:
  int rows, cols;
  unsigned char *data;

public:
  PGMImage() : rows(0), cols(0), data(0) {}

  ~PGMImage() { this->clear(); }

  void clear() {
    if (data) {
      delete[] data;
    }
    data = 0;
    rows = cols = 0;
  }

  int getRows() const { return (this->rows); }

  int getCols() const { return (this->cols); }

  void readFromFile(const string &fname) {
    char magic[2];
    ifstream ifile;
    int max;

    ifile.open(fname.c_str());
    if (!ifile) {
      throw(Exception("Could not open image file for reading"));
    }

    ifile.read(magic, 2 * sizeof(unsigned char));

    if (magic[0] != 'P' || magic[1] != '5') {
      throw(Exception("Image file not a PGM"));
    }

    ifile >> this->cols >> this->rows >> max;

    this->data = new unsigned char[this->rows * this->cols];

    ifile.read((char *)this->data,
               this->rows * this->cols * sizeof(unsigned char));

    ifile.close();
  }

  unsigned char operator()(const unsigned &r, const unsigned &c) const {
    return (this->data[r * this->cols + c]);
  }

  void writeToFile(const string &fname) const {
    ofstream ofile;

    ofile.open(fname.c_str());
    if (!ofile) {
      throw Exception("Could not open image for writing");
    }

    ofile << "P5\n" << this->cols << " " << this->rows << "\n255\n";

    ofile.write((char *)this->data,
                this->rows * this->cols * sizeof(unsigned char));

    ofile.close();
  }
};

class RadonSpace {
private:
  int numAngs;
  int numDists;
  double minDist;
  double maxDist;
  int *accumulator;

public:
  RadonSpace(const int &numA, const int &numD, const double &minD,
             const double &maxD)
      : numAngs(numA), numDists(numD), minDist(minD), maxDist(maxD),
        accumulator(0) {
    this->accumulator = new int[numA * numD];
    memset(this->accumulator, 0, numA * numD * sizeof(int));
  }

  ~RadonSpace() { this->clear(); }

  void clear() {
    if (accumulator) {
      delete[] accumulator;
    }
    accumulator = 0;

    numAngs = numDists = 0;
    minDist = maxDist = 0;
  }

  void computeSpace(const PGMImage &img) {
    unsigned rind;
    double cang, cdist, cs, sn;

    for (int a = 0; a < this->numAngs; a++) {
      cang = double(a) / (this->numAngs - 1) * M_PI;
      cs = cos(cang);
      sn = sin(cang);

      for (int d = 0; d < this->numDists; d++) {
        cdist =
            double(d) / (this->numDists - 1) * (this->maxDist - this->minDist) +
            this->minDist;

        for (unsigned cind = 0; cind < unsigned(img.getCols()); cind++) {
          rind = unsigned((cdist - cind * cs) / sn + 0.5);

          if (rind < unsigned(img.getRows())) {
            this->accumulator[a * this->numDists + d] += img(rind, cind);
          }
        }
      }
    }
  }

  void writeToFile(const string &fname) {
    ofstream ofile;
    int maxvotes;
    unsigned char nvotes;
    double scale;

    ofile.open(fname.c_str(), ios::binary);
    if (!ofile) {
      throw Exception("could not open file for writing");
    }

    ofile << "P5\n" << this->numDists << " " << this->numAngs << "\n255\n";

    maxvotes = accumulator[0];
    for (int i = 0; i < this->numAngs * this->numDists; i++) {
      if (maxvotes < this->accumulator[i]) {
        maxvotes = this->accumulator[i];
      }
    }

    scale = 255.0 / ((maxvotes <= 0) ? 1 : maxvotes);
    for (int a = 0; a < this->numAngs; a++) {
      for (int d = 0; d < this->numDists; d++) {
        nvotes =
            (unsigned char)(this->accumulator[a * this->numDists + d] * scale +
                            0.5);
        ofile.write((char *)(&nvotes), sizeof(unsigned char));
      }
    }

    ofile.close();
  }
};

int main(int argc, char **argv) {
  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <img.pgm>" << endl;
    return (1);
  }

  try {
    PGMImage img;
    img.readFromFile(argv[1]);

    RadonSpace acc(180, 3021, -1510, 1510);
    acc.computeSpace(img);

    acc.writeToFile("radonSpace.pgm");
  } catch (Exception &e) {
    cerr << e.what() << endl;
    return (1);
  }

  return (0);
}
