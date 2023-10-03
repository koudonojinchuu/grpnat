#include <iostream>
#include <vector>
#include <ranges>
#include <algorithm>

using namespace std;
using namespace ranges;

class Arrow_G { // an arrow is defined by each of its columns being the image of one of the base groups. Permutation matrix.
  public:
    vector<vector<unsigned>> m_coeffs;
    Arrow_G(vector<vector<unsigned>> coeffs): m_coeffs (coeffs) {};
    Arrow_G(unsigned n)
    {
      m_coeffs = vector<vector<unsigned>> (n, vector<unsigned> (n, 0));
    }
    // in a computer matrices are stored with lines the subunits, instead of columns as in maths notation or in TI92p etc.
    // This way the access would be still done (#line, #column).
    // The meaning of a column (math-vector) is unchanged.
    // It's just that the column needs to be accessed across several lines stored as different elements of the super cpp-vector.
    void print()
    {
      if (0 == m_coeffs.size()) {
        throw "coeffs matrix cannot contain 0 lines";
      }
      for (unsigned i = 0; i < m_coeffs.size(); ++i) {
        for (unsigned j = 0; j < m_coeffs[0].size(); ++j) {
        cout << m_coeffs[i][j];
        cout << " ";
       }
       cout << endl;
      }
      cout << "---" << endl;
    }
    vector<unsigned> apply(vector<unsigned> x)
    // x as a column vector
    {
      // check dims
      if (x.size() != m_coeffs.size()) {
        throw "incompatible sizes between x and coeffs matrix";
      }
      if (0 == m_coeffs.size()) {
        throw "coeffs matrix cannot contain 0 lines";
      }
      if (m_coeffs.size() != m_coeffs[0].size()) {
        throw "coeffs matrix needs to be square";
      }
      // do matrix-mul
      vector<unsigned> res (x.size(), 0);
      unsigned line_count = 0;
      for (auto line: m_coeffs) {
        unsigned col_count = 0;
        for (auto elt: line) {
          res[line_count] += elt*x[col_count];
          ++col_count;
        }
        ++line_count;
      }
      // output res
      return res;
    }
    vector<unsigned> inverse(vector<unsigned> x) // Arrow_G is assumed iso, so we can apply its inverse.
    {
      // check dims
      if (x.size() != m_coeffs.size()) {
        throw "incompatible sizes between x and coeffs matrix";
      }
      if (0 == m_coeffs.size()) {
        throw "coeffs matrix cannot contain 0 lines";
      }
      if (m_coeffs.size() != m_coeffs[0].size()) {
        throw "coeffs matrix needs to be square";
      }
      // do matrix-mul
      vector<unsigned> res (x.size(), 0);
      for (unsigned i: iota_view(0u, m_coeffs.size())) {
        for (unsigned j: iota_view(0u, m_coeffs.size())) {
          //res[i] += m_coeffs[i][j]*x[j]; // direct
          res[i] += m_coeffs[j][i]*x[j]; // inverse
        }
      }
      // output res
      return res;
    }
};

class Trnat {
  public:
    vector<vector<unsigned>> m_connected_components;
    vector<unsigned> m_invalidated_indices;
    unsigned m_nextindex = 1;

    Trnat(unsigned n): m_connected_components(vector<vector<unsigned>> (n, vector<unsigned> (n, 0))) {}

    void addstrut(unsigned i1, unsigned j1, unsigned i2, unsigned j2)
    {
      unsigned & index1 = m_connected_components[i1][j1];
      unsigned & index2 = m_connected_components[i2][j2];
      if (index1 == 0 && index2 == 0) {
        index1 = m_nextindex;
        index2 = m_nextindex;
        ++m_nextindex;
      } else if (index1 == 0) {
        index1 = index2;
      } else if (index2 == 0) {
        index2 = index1;
      } else if (index1 == index2) {
        // do nothing
      } else if (index1 > index2) {
        // set all index1-valued components to index2-value.
        for (auto & line: m_connected_components) for (auto & element: line)
          if (element == index1) element = index2;
        m_invalidated_indices.push_back(index1);
      } else if (index1 < index2) {
        // set all index2-valued components to index1-value.
        for (auto & line: m_connected_components) for (auto & element: line)
          if (element == index2) element = index1;
        m_invalidated_indices.push_back(index2);
      }
    }
};

void print(vector<unsigned> a)
{
  for (auto x: a) {
    cout << x << endl;
  }
  cout << "---" <<  endl;
}

void print(vector<vector<unsigned>> t)
{
  if (0 == t.size()) {
    throw "coeffs matrix cannot contain 0 lines";
  }
  for (unsigned i = 0; i < t.size(); ++i) {
    for (unsigned j = 0; j < t[0].size(); ++j) {
    cout << t[i][j];
    cout << " ";
   }
   cout << endl;
  }
  cout << "---" << endl;
}

// find all permutations, except identity, for a size n
// each given in the form of a vector <a b  c ...> meaning 0 turns into a, 1 into b, 2 into c, etc.
//vector<vector<unsigned>> find_permutations_except_id(unsigned n)
//{
//  vector<unsigned> used_array_bool (n, 0);
//
//}

vector<vector<unsigned>> find_perm(unsigned n)
{
  if (n == 0) {
    return {{}};
  //} else if (n == 1) {
  //  return {{0}};
  } else {
    vector<vector<unsigned>> result;
    auto vec_smaller_perms = find_perm(n - 1);
    for (auto & smaller_perm: vec_smaller_perms) {
      for (auto pos: iota_view (0u, n)) { // required to enforce "unsigned" by 0u or else template inference fails
        result.push_back(smaller_perm);
        auto & last = result.back();
        last.insert(last.begin()+pos, n - 1);
      }
    }
    return result;
  }
}

// convert a permutation into an arrow matrix
Arrow_G permut_to_arrowmat(vector<unsigned> p)
{
  unsigned n = p.size();
  Arrow_G result (n); // filled with 0
  for (auto i: iota_view(0u, n)) {
    for (auto j: iota_view(0u, n)) {
      if (p[i] == j) {
        result.m_coeffs[j][i] = 1;
      }
    }
  }
  return result;
}

// find all arrows (of dim n x n) which are permutations (<=> are bijections) and are not the identity
vector<Arrow_G> find_all_h(unsigned n)
{
  auto permutations = find_perm(n);
  // remove the identity (last element)
  permutations.erase(--permutations.end());
  vector<Arrow_G> result;
  for (auto x: permutations) result.push_back(permut_to_arrowmat(x));
  return result;
}

unsigned vec_to_scal(vector<unsigned> v) // only for vectors containing almost all zeros and one one.
{
  for (unsigned i: iota_view(0u, v.size()))
    if (v[i] == 1) return i;
  // error
  throw "this vector must be almost all zeros but one one";
}

vector<unsigned> scal_to_vec(unsigned n, unsigned i) // only for vectors containing almost all zeros and one one.
{
  vector<unsigned> res (n, 0); 
  res[i] = 1;
  return res;
}

int main()
{
  Arrow_G h1 ({{0,0,1},{0,1,0},{1,0,0}});
  h1.print();
  vector<unsigned> res = h1.apply({1,0,0});
  print(res);

  vector<vector<unsigned>> res2 = find_perm(2);

  for (auto x: res2) {
    for (auto y: x) {
      cout << y << ",";
    }
    cout << ";";
  }

  cout << endl;

  //for (auto x: iota_view(1,10)) cout << x << "#";
  auto res3 = find_perm(5);
  for (auto x: res3) {
   for (auto y: x) cout << y << ",";
   cout << ";";
  }

  cout << endl << endl;


  //for (auto h: allh) h.print();

  unsigned n = 10;

  Trnat tnat (n);
  auto allh = find_all_h(n);

  for (auto h: allh) {
    for (unsigned i: iota_view(0u, n))
      for (unsigned j: iota_view(0u, n)) {
        unsigned h_i = vec_to_scal(h.apply(scal_to_vec(n, i)));
        unsigned hr_j = vec_to_scal(h.inverse(scal_to_vec(n, j)));
        // tnat(h_i, j) = tnat(i, hr_j)
        tnat.addstrut(h_i, j, i, hr_j);
      }
  }

  print(tnat.m_connected_components);  


  return 0;
}
