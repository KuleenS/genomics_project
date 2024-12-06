/*
 * Four Russian class
 */

#include "FourRussians.h"
#include <omp.h> // Include OpenMP header for parallelization

FourRussians::FourRussians(string X, string Y, int tValue) {
  this->tValue = tValue;
  this->hash = Hash(tValue);
  this->stringA = X;
  this->stringB = Y;

  int numBlocks = getSizeBlocks();
  downOffsets = new uint8_t[numBlocks];
  rightOffsets = new uint8_t[numBlocks];
  blocks = new Block[numBlocks];

  editStrings();
}

void FourRussians::generate(int i, int tValue, vector<int>& counters) {
  if (i == tValue * 4) {
    int* topO = new int[tValue];
    int* leftO = new int[tValue];
    signed char* topS = new signed char[tValue];
    signed char* leftS = new signed char[tValue];

    string s[4];
    for (int k = 0; k < 4; ++k) {
      s[k] = "";
      for (int i = tValue * k; i < tValue * (k + 1); ++i) {
        if (k < 2) {
          s[k] += to_string(O[counters[i]]);
        } else {
          s[k] += E[counters[i]];
        }
      }
    }
    int hashed = hash.toHash(s[0], s[1], s[2], s[3]);

    for (int i = 0; i < tValue; i++) topO[i] = O[counters[i]];
    for (int i = tValue; i < tValue * 2; i++)
      leftO[i - tValue] = O[counters[i]];
    for (int i = tValue * 2; i < tValue * 3; i++)
      topS[i - tValue * 2] = E[counters[i]];
    for (int i = tValue * 3; i < tValue * 4; i++)
      leftS[i - tValue * 3] = E[counters[i]];

    Block b = Block(tValue, topS, leftS, topO, leftO);


    //#pragma omp critical
    {
      blocks[hashed] = b;
      downOffsets[hashed] = b.calcDownOffsets(tValue);
      rightOffsets[hashed] = b.calcRightOffsets(tValue);
    }

    delete[] topO;
    delete[] leftO;
    delete[] topS;
    delete[] leftS;

    return;
  } else {
    bool generate_letters = i >= tValue * 2;
    int k = generate_letters ? 4 : 3;

    for (int l = 0; l < k; ++l) {
      counters[i] = l;
      generate(i + 1, tValue, counters);
    }
  }
}

int FourRussians::calculateTValue(int lenA, int lenB) noexcept(true) {
  int len = max(lenA, lenB);
  double t = (log(len) / log(12)) / 2;

  return ((int)t + 1);
}

int** FourRussians::calculateEditMatrix() {
  int lenA = (int)stringA.size();
  int lenB = (int)stringB.size();

  getsubArrays();

  int** matrix = new int*[lenB / tValue];
  int lenb = lenB / tValue;
  int lena = lenA / tValue;

 //#pragma omp parallel for
  for (int h = 0; h < lenb; h++) {
    matrix[h] = new int[lena];
    for (int j = 0; j < lena; j++) {
      matrix[h][j] = 0;
    }
  }

  string sL, sT;
  for (int k = 0; k < tValue; k++) {
    sT = sT + "1";
    sL = sL + "1";
  }

  int leftS = substringB[0];
  int topS = substringA[0];

  matrix[0][0] = hash.toHash(sT, sL, topS, leftS);
  for (int k = 1; k < lena; k++) {
    topS = substringA[k];
    matrix[0][k] = hash.toHash(sT, rightOffsets[matrix[0][k - 1]], topS, leftS);
  }
  int initialtopS = substringA[0];

  // parallel for
  for (int i = 1; i < lenb; i++) {
    leftS = substringB[i];

    matrix[i][0] = hash.toHash(downOffsets[matrix[i - 1][0]], sL, initialtopS, leftS);

    for (int j = 1; j < lena; j++) {
      matrix[i][j] = hash.mergeHashes(downOffsets[matrix[i - 1][j]],
                                      rightOffsets[matrix[i][j - 1]],
                                      substringA[j], leftS);
    }
  }
  return matrix;
}

size_t FourRussians::getSizeBlocks() {
  size_t val = pow(3, tValue * 2) * pow(4, tValue * 2);
  return val;
}

void FourRussians::editStrings() {
  while (stringA.size() % tValue != 0)
    stringA = stringA.substr(0, stringA.size() - 1);
  while (stringB.size() % tValue != 0)
    stringB = stringB.substr(0, stringB.size() - 1);
}

int FourRussians::calculateMinDistance(int** matrix) {
  int min = 0;

  for (int j = 0; j < stringA.size() / tValue; j++) {
    min = min + blocks[matrix[stringB.size() / tValue - 1][j]].sumDown;
  }
  min = min + (int)stringB.size();
  return min;
}

void FourRussians::getsubArrays() {
  substringA = new int[stringA.size() / tValue];
  substringB = new int[stringB.size() / tValue];
  for (int i = 0; i < stringA.size() / tValue; i++) {
    substringA[i] = 0;
    for (int j = 0; j < tValue; ++j) {
      substringA[i] =
          substringA[i] * 4 + hash.letterToNum(stringA[i * tValue + j]);
    }
  }
  for (int i = 0; i < stringB.size() / tValue; i++) {
    substringB[i] = 0;
    for (int j = 0; j < tValue; ++j) {
      substringB[i] =
          substringB[i] * 4 + hash.letterToNum(stringB[i * tValue + j]);
    }
  }
}

void FourRussians::parallelPreProcessing() {
  E[0] = 'A';
  E[1] = 'C';
  E[2] = 'T';
  E[3] = 'G';
  O[0] = -1;
  O[1] = 0;
  O[2] = 1;

  vector<thread> threads;

  for (int c0 = 0; c0 < 3; ++c0) {
    for (int c1 = 0; c1 < 3; ++c1) {
      for (int c2 = 0; c2 < 3; ++c2) {
        auto code = [this, c0, c1, c2]() {
          vector<int> comb;
          comb.resize(FourRussians::tValue * 4);
          comb[0] = c0;
          comb[1] = c1;
          comb[2] = c2;
          generate(3, tValue, comb);
        };
        threads.push_back(thread(code));
      }
    }
  }

  for (auto& t : threads) {
    t.join();
  }
}
