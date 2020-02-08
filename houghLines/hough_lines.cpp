#include <exception>
#include <fstream>
#include <iostream>
#include <math.h>
#include <string.h>
#include <vector>

using namespace std;

const int ANGULAR_OVERSAMP = 100;

class Exception : public exception {
private:
  const char *message;

public:
  Exception(const char *m) : message(m) {}

  virtual const char *what() const throw() { return (this->message); }
};

class Point2 {
public:
  int r, c;
};

class HoughSpace {
private:
  int numAngs;
  int numDists;
  double minDist;
  double maxDist;
  int *accumulator;

public:
  HoughSpace(const int &numA, const int &numD, const double &minD,
             const double &maxD)
      : numAngs(numA), numDists(numD), minDist(minD), maxDist(maxD),
        accumulator(0) {
    this->accumulator = new int[this->numAngs * this->numDists];
    memset(this->accumulator, 0, this->numAngs * this->numDists * sizeof(int));
  }

  ~HoughSpace() { this->clear(); }

  void clear() {
    if (accumulator) {
      delete[] accumulator;
    }
    accumulator = 0;

    numAngs = numDists = 0;
    minDist = maxDist = 0;
  }

  void addVote(const int &r, const int &c) {
    double cang, cdist;
    int aind, dind;

    for (int a = 0; a < this->numAngs * ANGULAR_OVERSAMP; a++) {
      cang = double(a) / (this->numAngs * ANGULAR_OVERSAMP - 1) * M_PI;
      cdist = c * cos(cang) + r * sin(cang);

      aind = int(double(a) / (this->numAngs * ANGULAR_OVERSAMP - 1) *
                     (this->numAngs - 1) +
                 0.5);

      dind = int((cdist - this->minDist) / (this->maxDist - this->minDist) *
                     (this->numDists - 1) +
                 0.5);

      this->accumulator[aind * this->numDists + dind]++;
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
  ifstream ifile;
  int row, col;
  vector<Point2> pts;
  Point2 pt;

  if (argc != 2) {
    cerr << "usage: " << argv[0] << " <pts.txt>" << endl;
    return (1);
  }

  ifile.open(argv[1]);
  if (!ifile) {
    cerr << "Error: could not open points file for reading" << endl;
    return (1);
  }

  while (ifile >> row >> col) {
    pt.r = row;
    pt.c = col;

    pts.push_back(pt);
  }

  ifile.close();

  try {
    HoughSpace acc(180, 285, -142, 142);
    for (vector<Point2>::iterator it = pts.begin(); it != pts.end(); it++) {
      acc.addVote(it->r, it->c);
    }

    acc.writeToFile("houghSpace.pgm");
  } catch (Exception &e) {
    cerr << e.what() << endl;
    return (1);
  }

  return (0);
}
