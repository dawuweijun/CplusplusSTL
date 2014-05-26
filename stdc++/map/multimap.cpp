// erasing from multimap
#include <iostream>
#include <map>

void print(std::multimap<char, int>& mmap)
{
    std::cout << "Now all the content are:\n";
    std::multimap<char,int>::iterator it; 
  for (it=mmap.begin(); it!=mmap.end(); ++it)
      std::cout << (*it).first << " => " << (*it).second << '\n';
}
int main ()
{
  std::multimap<char,int> mymultimap;

    // insert some values:
  mymultimap.insert(std::pair<char,int>('a',10));
  mymultimap.insert(std::pair<char,int>('b',20));
  mymultimap.insert(std::pair<char,int>('b',30));
  mymultimap.insert(std::pair<char,int>('c',40));
  mymultimap.insert(std::pair<char,int>('d',50));
  mymultimap.insert(std::pair<char,int>('d',60));
  mymultimap.insert(std::pair<char,int>('e',70));
  mymultimap.insert(std::pair<char,int>('f',80));

  std::multimap<char,int>::iterator it = mymultimap.find('b');
  print(mymultimap);

  mymultimap.erase (it);                     // erasing by iterator (1 element)
  print(mymultimap);

  mymultimap.erase ('d');                    // erasing by key (2 elements)
  print(mymultimap);

  it=mymultimap.find ('e');
  mymultimap.erase ( it, mymultimap.end() ); // erasing by range
  print(mymultimap);
  return 0;
}
