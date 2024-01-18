/*
  SWARM

  Copyright (C) 2012-2024 Torbjorn Rognes and Frederic Mahe

  This program is free software: you can redistribute it and/or modify
  it under the terms of the GNU Affero General Public License as
  published by the Free Software Foundation, either version 3 of the
  License, or (at your option) any later version.

  This program is distributed in the hope that it will be useful,
  but WITHOUT ANY WARRANTY; without even the implied warranty of
  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
  GNU Affero General Public License for more details.

  You should have received a copy of the GNU Affero General Public License
  along with this program.  If not, see <http://www.gnu.org/licenses/>.

  Contact: Torbjorn Rognes <torognes@ifi.uio.no>,
  Department of Informatics, University of Oslo,
  PO Box 1080 Blindern, NO-0316 Oslo, Norway
*/

#include <cassert>
#include <string>
#include <vector>


auto compress_alignment_to_cigar(std::vector<char> const & input,
                                 std::string & destination) -> void {
  assert(not input.empty());

  // run-length encoding (RLE) algorithm
  auto is_first = true;
  auto previous_char = '\0';
  auto count = 0UL;
  for (const auto current_char: input) {
    if (is_first) {
      previous_char = current_char;
      ++count;
      is_first = false;
      continue;
    }
    if (current_char == previous_char) {
      ++count;
      continue;
    }
    if (count > 1) {
      destination.append(std::to_string(count));
      count = 1;
    }
    destination.push_back(previous_char);
    previous_char = current_char;
  }
  
  // last item
  if (count > 1) {
    destination.append(std::to_string(count));
  }
  destination.push_back(input.back());
}


// tests already covered in swarm-tests:

//  - singleton: M -> M
//  - monotonous: MMM -> 3M
//  - initial singleton: DMM -> D2M
//  - terminal singleton: MMD -> 2MD
//  - double-digit count: MMMMMMMMMM -> 10M
//  - triple-digit count: MMM...MMMM -> 100M
//  - quadruple-digit count: MMM...MMMM -> 1000M


// tests (compiler explorer):

// int main() {
//     std::string cigar;

//     // alternated: MIMIMI -> MIMIMI
//     std::vector<char> alternated{'M', 'I', 'M', 'I', 'M', 'I'};
//     compress_alignment_to_cigar(alternated, cigar);
//     assert(cigar == "MIMIMI");
//     std::cout << cigar << '\n';
//     cigar.clear();

//     // double alternated: MMIIMMIIMMII -> 2M2I2M2I2M2I
//     std::vector<char> double_alternated{'M', 'M', 'I', 'I', 'M', 'M',
//                                         'I', 'I', 'M', 'M', 'I', 'I'};
//     compress_alignment_to_cigar(double_alternated, cigar);
//     assert(cigar == "2M2I2M2I2M2I");
//     std::cout << cigar << '\n';
//     cigar.clear();
// }

// see:
// https://godbolt.org/#z:OYLghAFBqd5QCxAYwPYBMCmBRdBLAF1QCcAaPECAMzwBtMA7AQwFtMQByARg9KtQYEAysib0QXACx8BBAKoBnTAAUAHpwAMvAFYTStJg1DIApACYAQuYukl9ZATwDKjdAGFUtAK4sGIAMxcpK4AMngMmAByPgBGmMQgAKyJpAAOqAqETgwe3r4BQemZjgJhEdEscQnJtpj2JQxCBEzEBLk%2BfoG19dlNLQRlUbHxSSkKza3t%2BV3j/YMVVaMAlLaoXsTI7Bzm/uHI3lgA1Cb%2BbqIKSq0n2CYaAII7eweYx6dO48SYrNe3D2a7DH2XiOJzcH3CwB%2B90egOerzcADdMA4SFCHvcmF4iIc0CxUp8LgB9MR4YAMNiCQlEQnIUktCDjdAgEBIlHEUHIBAta44gTjcwANkO4VSWNIv0Okql0plsrlMsZzPBRkFhyw43CTAaS0OAFoeQjUHh0McAOxWe7SvBUQ4QEVYgB0mDxBAAnhAljqTOaJXLPgR1gwTha7tLvQARX6%2ByUAehjh2IXgYuvoRgICEOjDQ%2BCMhzEwBIhAQLAU0bzWNQwoUhJoxHGr3DhwIicwwbLmOx%2BMwCKcXmrnJaDcOYG2iTcGhHbctUo7lbQSYIQ40chCU9DUv4xFtaAY9dnOPWn0pA63IGFDFFBC9Punsuttrw1dr42vIfl0q7PbW/a5W5OjeQQ9GAIGlfzXd8pWsax50EcCIKrGs8DrRd/0OKgxCUOCIJ3RwGC8Vt/DfOUIzLK0bQgQDiCPECTwbVDP17H8WlfUjiMsaC1lgwjWNlHDwnwrCZRI28ZXvCjOJQ/xsEOLgWJE%2BV1VwrVsgdJhUlSVwGQIJkQGpZVgHEhdPUE%2BUYMkxsuBMyDTUjeSpUUzUGgdUUFAQQkYiYZAAGsIAY79QOYqzJT8vsAr/fwAKA48wO4kThPRddY3jAx60IZ0yzEszXmk2SzSI2UHOYJy1I0hh0C0nS9ObCFDMEYzYsSs1bMawrlIEZy%2BzcjzvLtC9HW6nz6pDeLfnCRcWCYcIPTystFRAfScTpdkGrLOM81oAh4iKzAmUOABZABJQ7Dr1HljqOg7Zu05lWSIZazhi6SxE24htvQb0LBHPaR1IYcOAOn6/u%2B3g/oBkGvsBkcwZIhrpVxLsiRJMkKRA6laWAelnq2rUdt%2B9HAthmcLniAgKKWuiIuOMwzHOw7zDMJZBLm%2BdJLcUFFox8K2dOP6TDHBhJ0JyV8eIB19i%2BYgPTXVb43QNYYnodaXres89sOo71bVg6DtOqTDhpswDoNo29sNq6dNu1FThPHk5a8BXMGJDbsc297zQh8GOGB36och/6/e9oGfp4%2BDQ4g33PbBn2vYDv2o9BycbME%2BGCWrJHyWAqlUBpJaIDth2neVnH0DxpbGaFvNidaMnOYpxt6eNxvTcu6ny/yw5mbWVn2ZF%2BF2cnfnBfbkWxfoek26jGyOBWWhOESXg/A4LRSFQTg2fYywO7WDYXh2HhSAITRp5WLykg0fROEkBej5XzheAUEBz8Ppfp9IOBYBgRAUFQPE6HichKDwz/gkYAChmCpFcqgAgfA6AvQfhAGIN8YiamIK6Tg%2B9kHMFQQAeRiNoZEz9l64hRtghgtA0Ev1IFgGIXhgBuDELQB%2B3BeBYAmkYcQlD8CfAcHgJETDl6YFUMiLEWxl5jTqDfWgeAYjEBaK6DwWAb7VRYOg1%2BVADCgIAGp4EwAAd2waVVRMhBAiDEOwKQxj5BKDUDfXQQQDBGBQFBSw%2BhpEP0gCsVAqQGhMN1Iyf8pgN4WC4FwPU2D/C8FQEiKixpMDuI9N0Ah2QXBlSmH4f4wQyrzGGNUc%2BRQsgCDSSADJ%2BSGjZMqCMRI587BJIEH0SYngOjFP8IknhdSJgDHCEMCpuTbAdKKRk2YrRymLCqSsBQ29NgSBnnPa%2BlDV4cEOKoAAHAKXUApJCHGAMgZAtpmxJi8jqCAuBCAkCpoEJYvBCHH1IKfKpF8OBX1ICo%2B5i9l4LPvo/A%2BR8Vjvy/l3S8ADxK/3oMQSIrAtgrLWRsrZOzDiSAdMsh089gj4DusaPQ/ATGiHEBYrFViVDqEoXY0gujZGpFUTMjgKK3mRM4NgrEl5DioBtFC9Zmztm7IgPshghzbQeBBfEc5skrk/JWAgL4WAEgJNno83gLzz60tvhwT5T8xW3LPg8iJpAlUfO%2BS/T0DyzBzPeXffVWhDXRMyM4SQQA%3D
