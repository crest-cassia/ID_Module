# -*- coding: utf-8 -*-
require 'json'
# RSRuby
require 'rsruby'
r = RSRuby.instance

# このプログラムは、R のプログラムをそのままの形で ruby に書き換えたものである。
# 評価関数 f() を、OACISを呼び出すことなく実行し、交互作用判定を実行する。
# 
# OACISで実行する場合は、f()の呼び出し順が R と異なる。
#  R では、置換直交表を、1行めから最後の行まで f()に適用して１回の結果を得て、
#  それを repeat 回繰り返す。
#  OACIS では、置換直交表の１行について repeat 回実行し、それを置換直交表の
#  すべての行について行う。
# このため、R と OACISでは、得られる計算結果の「計算順」が異なることになる。
#
# このプログラムには、OACISの計算順をシミュレートする関数を用意してある。
#  (デフォルトでコメントアウトとしてある)
# これらの関数は、OACISから得られた計算結果を処理する関数をデバッグする目的で
# 使用したが、OACISの計算順をシミュレートすることにも使用可能である。
#
# このプログラムは、判定結果を標準出力に表示する(ファイルのは出力しない)。
#
# OACIS版と違って、f() は、このプログラム中に定義する。
# 水準値は _level_value.json ファイルから読み出す。 
# "open" で検索すると、jsonファイルを読んでいる箇所がわかるはず。


class Array 
  def average
    inject(0.0) {|sum, i| sum += i }/size
  end

  def variance
    ave = average
    inject(0.0) {|sum, i| sum += (i - ave)**2 }/(size-1)
  end

  def standard_devitation
    Math::sqrt(variance)
  end
end


def norm_rand(mu=0,sigma=1.0)
  r1,r2=Random.rand(),Random.rand();
  Math::sqrt(-1*Math::log(r1))*Math::cos(2*Math::PI*r2)*sigma+mu
end


def f(*v)
#   return 10*v[0] + 10*v[1] + 10*v[2]*v[3]

     #Example 0
#     10*v[0]*v[1]

     #Example 1, 2
     #100*v[0]+100*v[1]+100*v[3]*v[4]+norm_rand(mu=0,sigma=0.3*v[5]+1)

     # Example 5
#    1000000*v[0]+1000000*v[1]+1000000*v[3]*v[4]+norm_rand(mu=0,sigma=1.0)

     # Example 3
#     100*v[0]*100*v[1]+100*v[3]*v[4]+norm_rand(mu=0,sigma=0.3*v[5]+1)+100*v[10]+100*v[0]*v[12]**2*v[13]+v[8]*v[9]*v[10]*v[11]+100*v[2]+0.01*v[14]**2;

     # Example 4
#     100*v[63]*100*v[117]+100*v[115]*100*v[64]+100*v[17]*100*v[29]+100*v[16]*100*v[89]+100*v[47]*100*v[84]*100*v[21]+100*v[92]*100*v[54]+100*v[114]*100*v[22]+100*v[15]*100*v[126]+100*v[62]*100*v[69]+100*v[126]*100*v[26]*100*v[83]+100*v[35]*100*v[12]+100*v[122]*100*v[34]*100*v[62]+100*v[95]*100*v[67]+100*v[124]*100*v[121]*100*v[39]+100*v[50]*100*v[23]+100*v[58]*100*v[44]+100*v[126]*100*v[10]+100*v[45]*100*v[57]+100*v[84]*100*v[12]+100*v[109]*100*v[29]*100*v[33]+100*v[19]*100*v[97]+100*v[11]*100*v[65]*100*v[115]*100*v[31]+100*v[24]*100*v[4]+100*v[45]*100*v[103]+100*v[20]*100*v[101]+100*v[89]*100*v[103]*100*v[111]+100*v[43]*100*v[1]+100*v[36]*100*v[21]+100*v[20]*100*v[60]+100*v[74]*100*v[22]*100*v[63]*100*v[78]+100*v[13]*100*v[21]*100*v[95]+100*v[108]*100*v[3]+100*v[116]*100*v[77]+100*v[20]*100*v[69]+100*v[36]*100*v[11]+100*v[91]*100*v[75]*100*v[105]*100*v[82]+100*v[46]*100*v[54]+100*v[92]*100*v[0]*100*v[65]+100*v[21]*100*v[63]*100*v[98]+100*v[43]*100*v[65]*100*v[66]+100*v[107]*100*v[31]*100*v[115]+100*v[12]*100*v[33]*100*v[72]+100*v[31]*100*v[105]*100*v[22]+100*v[75]*100*v[47]+100*v[2]*100*v[83]*100*v[71]+100*v[93]*100*v[92]+100*v[76]*100*v[78]+100*v[106]*100*v[54]+100*v[119]*100*v[9]+100*v[105]*100*v[8]*100*v[35]+norm_rand(mu=0,sigma=1.0)

     # Example 4 ver2
     100*v[30]*100*v[50]+100*v[26]*100*v[9]+100*v[4]*100*v[50]+100*v[35]*100*v[19]+100*v[62]*100*v[56]+100*v[53]*100*v[51]+100*v[29]*100*v[60]+100*v[58]*100*v[60]+100*v[7]*100*v[49]+100*v[10]*100*v[49]+100*v[49]*100*v[18]+100*v[40]*100*v[27]+100*v[0]*100*v[35]+100*v[47]*100*v[1]+100*v[34]*100*v[46]+norm_rand(mu=0,sigma=1.0)


     # Example 4 ver3
#     100*v[7]*100*v[26]+100*v[4]*100*v[15]+100*v[18]*100*v[1]+100*v[4]*100*v[20]+100*v[11]*100*v[15]+100*v[0]*100*v[1]+100*v[4]*100*v[27]+100*v[7]*100*v[1]+100*v[21]*100*v[4]+100*v[14]*100*v[2]+100*v[24]*100*v[13]+100*v[7]*100*v[22]+100*v[24]*100*v[9]+100*v[4]*100*v[16]+100*v[23]*100*v[29]+norm_rand(mu=0,sigma=1.0)


end

### Function definitions.

def rev_orthogonal_table(ortable, check)
   wk = Array.new(ortable.length)
   (wk.length).times {|i|
      wk[i] = ortable[i][check]
   }
#   print "wk = #{wk}\n"
   lev = wk.uniq.sort

   rev_orth_list = Array.new(lev.length, 0)
   (lev.length).times {|j|
      qqq = Marshal.load(Marshal.dump(ortable))
      (ortable.size).times {|i|
         qqq[i][check] = lev[j]
#        p "qqq[i] = #{qqq[i]}"
      }
      rev_orth_list[j] = qqq
   }
   rev_orth_list
end


### These routines made for OACIS module operatons. ###
def or_make_data_set(rev_orth1, v2)
   param = Marshal.load(Marshal.dump(rev_orth1))
#   print "param"
#   p param
#   print "v2"
#   p v2
   (param.length).times {|k|
      (param[k].length).times {|j|
         (param[k][j].length).times {|i|
#            print "v2[ #{k} ] = #{v2[k]}\n"
#            print "param[#{k}][#{j}][#{i}] = #{param[k][j][i]}\n"
            param[k][j][i] = v2[i][param[k][j][i].to_i - 1]
         }
      }
   }
#   p param
   param
end


def or_get_parameter_set(v2=vv2, nfactors=1, plan=nil)
   rev_orth = Array.new(nfactors, 0)
   param_set = Array.new(nfactors, 0)
   nfactors.times {|factor1|
#      data1[factor1] = make_data_set(v1, v2, factor1, plan)
      rev_orth[factor1] = rev_orthogonal_table(plan, factor1)
      param_set[factor1] = or_make_data_set(rev_orth[factor1], v2)
   }
#   print "rev_orth\n"
#   p rev_orth
#   print "param_set\n"
#   p param_set
   return rev_orth, param_set
end

def dump_array(array, size, title)
   size.times {|k|
      (array[k].length).times {|j|
         print "#{title}[#{k}][#{j}]\n"
         (array[k][j].length).times {|i|
            print "#{array[k][j][i]}\n"
         }
      }
   }
end


def count_param(array, size)
   count = 0
   size.times {|k|
      (array[k].length).times {|j|
         (array[k][j].length).times {|i|
            count += 1
         }
      }
   }
   count
end


def or_use_orthogonal_table(fun, params, repeat)
   number_of_experiment = params.length
   print "number_of_experiment = #{number_of_experiment}\n"
   out = Array.new(number_of_experiment, 0).map{Array.new(repeat,0)}
   number_of_experiment.times {|i|
      print "params = #{params[i]}\n"
      repeat.times {|j|
         out[i][j] = fun.call(*params[i])
      }
#print "out[#{i}] = #{out[i]}\n"
   }
   out
end


def or_sampling_data(fun, param_set, repeat=1)
   factors = param_set.length
#   output = Array.new(factors, 0)
   output = []

#   print "factors= #{factors}\n"
   factors.times {|j|
      levels = param_set[j].length
#      print "levels= #{levels}\n"
#      out = Array.new(levels, 0)
      levels.times {|i|
         x2 = param_set[j][i]
#         print "x2= #{x2}\n"
#         out[i] = use_orthogonal_table(fun, x2, repeat)
         result = or_use_orthogonal_table(fun, x2, repeat)
         output << result
      }
#      output[j] = out
   }
#print "output #{output}\n" 
   output
end

def arrange_data_format(nfactors, nlevels, rev_orth, data1, repeat)
   data_set = Array.new(nfactors).map{Array.new(2, 0)}
   num_experiment = rev_orth[0][0].length
   print "nfactors = #{nfactors}\n"
   print "nlevels = #{nlevels}\n"
   print "repeat = #{repeat}\n"
   print "num_experiment = #{num_experiment}\n"

# print "data1 #{data1}\n"
   counter = 0
   nfactors.times {|factor|
      result2 = Array.new(nlevels[factor],0)
      (nlevels[factor]).times {|level|
         result1 = Array.new(repeat).map{Array.new(num_experiment,0)}
         repeat.times {|rep|
            num_experiment.times {|nexp|
               result1[rep][nexp] = data1[counter][nexp][rep]
            }
            result2[level] = result1
         }
         counter += 1
         data_set[factor][0] = result2
      }
      data_set[factor][1] = rev_orth[factor]
   }
   data_set
end

def conv_flat_data(data, nexp, nlevels)
   tlevel = nlevels.inject(:+)
   counter = 0
   data1 = Array.new(tlevel)
   tlevel.times {|cnt|
      data1[cnt] = []
      nexp.times {|loop|
         data1[cnt] << data[counter + loop]
      }
      counter += nexp
   }
   data1
end
# End of OACIS interface functions.


## Debug function definitions.
## Calculate same order as R program. 
def sampling_data(fun, v1=vv1, v2=vv2, repeat1=30, plan=NULL)
  len = v1.length
  data1 = Array.new(len, 0)
  len.times {|factor1|
    data1[factor1] = make_data_set(fun, v1, v2, factor1, repeat1, plan)
  }

  print "sampling_data data is \n"
  p data1

  data1
end


def make_data_set(fun, v1, v2, check=1, repeat1=1, plan=NULL)
  plan_base = Marshal.load(Marshal.dump(plan))
  rev_plan = rev_orthogonal_table(plan_base,check)

  print "make_data_set: rev_plan\n"
  p rev_plan

  q_rev = Array.new(rev_plan.length, 0)
  (rev_plan.length).times {|i|
    q_rev[i] = repeat_experiment_b(fun, v2, rev_plan[i], repeat1)
  }

  list1 = [q_rev, rev_plan]

#  print "make_data_set: list1\n"
#  p list1
  list1
end


def repeat_experiment_b(fun, v2, orTable, repeat1=10)
  output = Array.new(repeat1, 0)

  repeat1.times {|i|
    qq1 = use_orthogonal_table(fun, v2, orTable, x=1.0)
    output[i] = qq1
  }
  output
end


def use_orthogonal_table(fun, v2, orTable, x1=1.0)
  number_of_experiment = -1
  #直交表の行数の取得	
  number_of_experiment = orTable.length

  printf "number_of_experiment = %d\n", number_of_experiment

  out = Array.new(number_of_experiment, 0)
  number_of_experiment.times {|i|
    x2 = orTable[i]
#    print "x2 = "
#    p x2

    v = Array.new(x2.length, 0)
    (x2.length).times {|j2|
      v[j2]=v2[j2][x2[j2].to_i - 1]
#      printf "j2 = %d, v2[%d][%d]\n" , j2, j2, x2[j2]-1
    }
#    print "v = "
#    p v
    out[i] = fun.call(*v)
  }
  p out
  out
end

## End of debug function definitions.



def dump_data_set(data_set, nfactors, nlevels, repeat)
   nfactors.times {|factor|
      print "factor = #{factor}\n"
      (nlevels[factor]).times {|level|
         print "  level = #{level}\n"
         repeat.times {|nrep|
            print "    repeat #{nrep}  "
            print "#{data_set[factor][0][level][nrep]}\n"
         }
         print "    Rev orth\n"
         (data_set[factor][1][level].length).times {|nexp|
            print "      #{data_set[factor][1][level][nexp]}\n"
         }
      }
   }
end


def fout_data_set(data_set, nfactors, nlevels, repeat)
  open("dataset.out", "w") { |f|
    nfactors.times {|factor|
      print "factor = #{factor}\n"
      (nlevels[factor]).times {|level|
        print "  level = #{level}\n"
        repeat.times {|nrep|
          print "    repeat #{nrep}  "
          print "#{data_set[factor][0][level][nrep]}\n"
        }
        print "    Rev orth\n"
        (data_set[factor][1][level].length).times {|nexp|
          print "      #{data_set[factor][1][level][nexp]}\n"
        }
      }
    }
  }
end


def fout_orth(plan)
  open("orth.out", "w") { |f|
    (plan.length).times {|i| print "#{plan[i]} \n" }
  }
end


def fout_revorth(rev_orth, nfactors)
  open("revorth.out", "w") { |f|
    nfactors.times {|k|
      (rev_orth[k].length).times {|j|
        print "[#{k}][#{j}]\n"
        (rev_orth[k][j].length).times {|i|
          print "#{rev_orth[k][j][i]}\n"
        }
      }
    }
  }
end


def finding_interaction_all_c(data_set)
  len = data_set.length
  list1 = Array.new(len, 0)
  len.times {|i|
    list1[i] = finding_interaction_c(data_set, i)
  }
#  print "finding_interaction_all_c returns \n"
#  p list1
  list1
end


def finding_interaction_c(data_set, check)

  #要因番号checkのデータを取り出す
  data_set_i = data_set[check]

  #前処理：実験の各行の繰り返し分散の計算 
  noise_data_i = diff_epsilon_c(data_set_i)

  #前処理:各要因をある水準に固定したときの誤差分散の実験行平均の計算
  error_i = average_epsilon_c(noise_data_i)

  #誤差の交互作用の検査
  error_test = test_interaction_epsilon_c(noise_data_i)

  #見せかけの効果の検査
  pseudo_test = test_pseudo_effect_c(data_set_i, epsilon=error_i)

  #要因の交互作用の検査
  confounding_test = test_effect_confounding_c(data_set_i, epsilon=error_i)

  #結果を返す
  results_i = [error_test, pseudo_test, confounding_test]
  results_i
end


def diff_epsilon_c(data_set_i)
  e_rev = Array.new(data_set_i.length, 0)
  d2 = data_set_i[0]
#  0.upto(data_set_i[0].length-1) {|i|
  (d2.length).times {|i|
    e_rev[i] = calc_epsilon_c(d2[i])
  }
  e_rev
end


def calc_epsilon_c(data_i_j)
  sd1 = Array.new(data_i_j[0].length, 0)
  (data_i_j[0].length).times {|k|
    q = Array.new(data_i_j.length, 0)
    (data_i_j.length).times {|repeat_no|
      q[repeat_no] = data_i_j[repeat_no][k]
    }
    sd1[k]=q.standard_devitation
  }
  sd1
end


def average_epsilon_c (diff_epsilon_b)
  q7 = diff_epsilon_b
#  ss = Array.new(diff_epsilon_b.length, 0)
  ss = Array.new(q7.length, 0)
  (q7.length).times {|i|
    ss[i] = q7[i].average
  }
  ss
end


def test_interaction_epsilon_c(data_set_epsilon=noise_data_i)
  d_rev = data_set_epsilon

  t = Array.new(d_rev.length).map{Array.new(d_rev.length,0)}
  (d_rev.length).times {|j|
    (d_rev.length).times {|i|
      sd1 = d_rev[j].average
      sd2 = d_rev[i].average

      print "sd1 sd2 #{sd1}, #{sd2}\n"
#p d_rev[j]
#p d_rev[i]

      threshold = 2
      if [sd1, sd2].min != 0 then
        if ( (sd1 - sd2).abs/[sd1, sd2].min > threshold) then
          t[i][j] = TRUE
print "answer0 TRUE: #{i+1} #{j+1}\n"
        else
          t[i][j] = FALSE
print "answer0 FALSE: #{i+1} #{j+1}\n"
        end
      else
        if sd1 == sd2 then
          t[i][j] = FALSE
print "answer0 FALSE: #{i+1} #{j+1}\n"
        else
          t[i][j] = TRUE
print "answer0 TRUE: #{i+1} #{j+1}\n"
        end
      end
    }
  }
  t
end

def test_pseudo_effect_c(data_set, epsilon)
  d_rev = data_set[0]
#print "d_rev = #{d_rev}\n"
#print "d_rev.len = #{d_rev.length}\n"
  t1 = Array.new(d_rev.length).map{Array.new(d_rev.length,0)}

# variable "v_max" is not used.

  (d_rev.length).times {|i|
    (d_rev.length).times {|j|
#      res = preprocessing(d_rev, i, j)
#      vv1 = res[0]
#      vv2 = res[1]
      vv1, vv2, dmy = preprocessing(d_rev, i, j)
#print "xxx vv1 = #{i}, #{j}, #{vv1}\n"
#print "xxx vv2 = #{i}, #{j}, #{vv2}\n"

      threshold = 2
      sd = Array.new(vv1.length,0)
      (sd.length).times {|k|
        sd[k] = (vv1[k] - vv2[k])**2
      }
#print "xxx sd = #{sd} \n"
#print "xxx average = #{sd.average} \n"
#print "sum = #{sd.inject{|sum, n| sum + n}} \n"
      tmp1 = Math::sqrt( sd.average )
      tmp2 = threshold * Math::sqrt((epsilon[i]**2 + epsilon[j]**2))
      print ("answer L R = #{tmp1} #{tmp2}\n")
      if Math::sqrt( sd.average ) >= threshold * Math::sqrt((epsilon[i]**2 + epsilon[j]**2)) then
          t1[i][j] = TRUE
print "answer TRUE: #{i+1} #{j+1}\n"
      else
          t1[i][j] = FALSE
print "answer FALSE: #{i+1} #{j+1}\n"
      end
    }
  }
  t1
end


def test_effect_confounding_c(data_set, epsilon)
  d_rev = data_set[0]

# variable "v_max" is not used.
  t2 = Array.new(d_rev.length).map{Array.new(d_rev.length,0)}
  (d_rev.length).times {|i|
    (d_rev.length).times {|j|
#      res = preprocessing(d_rev, i, j)
#      vv1 = res[0]
#      vv2 = res[1]
      vv1, vv2, dmy = preprocessing(d_rev, i, j)
      threshold = 2

      vv0 = Array.new(vv1.length, 0)
      (vv1.length).times {|k|
        vv0[k] = vv1[k] - vv2[k]
      }

      if ( vv0.standard_devitation >= threshold * Math::sqrt((epsilon[i]**2 + epsilon[j]**2))) then
        t2[i][j] = TRUE
#        print "answer TRUE: #{i}, #{j}\n"
      else
        t2[i][j] = FALSE
#        print "answer FALSE: #{i}, #{j}\n"
      end
    }
  }
  t2
end

def preprocessing(d_rev, i, j)
  q1 = d_rev[i]
  q2 = d_rev[j]

  vv1 = Array.new(q1[0].length, 0)
  vv2 = Array.new(q1[0].length, 0)

#  print "q1 start\n"
#  p q1
#  p q2
#  print "q1 end\n"
  (q1[0].length).times {|kk|
    v1 = Array.new(q1.length, 0)
    v2 = Array.new(q1.length, 0)
    (q1.length).times {|k|
#      v1[k] = q1[k][kk].flatten
#      v2[k] = q2[k][kk].flatten
      v1[k] = q1[k][kk]
      v2[k] = q2[k][kk]
    }
#    print "v1 #{v1}\n"
#    print "v2 #{v2}\n"
    vv1[kk] = v1.average
    vv2[kk] = v2.average
  }

  v3 = Array.new(vv1.length, 0)
  (vv1.length).times {|i|
    v3[i] = (vv1[i] - vv2[i]).abs
  }
#  print "vv1 #{vv1}\n"
#  print "vv2 #{vv2}\n"
#  print "answer i-j #{i} #{j}"
  return [vv1, vv2, v3.max]
end



def check_interaction_all_c(data_set_set)
  mle2 = Array.new(data_set_set.length)
  (mle2.length).times {|i|
    mle2[i] = check_interaction_c(data_set_set[i], i)
  }
  mle2
end

def check_interaction_c(data_set_i, check=1)
  table_list = data_set_i[1]
  pseudo_effect = Array.new(table_list.length, 0)

  (table_list.length).times {|level|
    print "level = #{level}\n"
#    table_factor_check_is_i = as_list_orthogonal_table(table_list[level].transpose)
#print "table_list[level] = #{table_list[level]}\n"
#    table_factor_check_is_i = table_list[level].transpose
    table_factor_check_is_i = table_list[level]
#print "table_factor_check_is_i = #{table_factor_check_is_i}\n"
    data_set_i_j = data_set_i[0][level]
    v = calc_effects_c_using_table(data_set_i_j, table_factor_check_is_i)
#print "v = #{v}\n"
    pseudo_effect[level] = v
  }

print "pseudo_effect = #{pseudo_effect}\n"

  mse_s = Array.new(pseudo_effect.length, 0)

  (pseudo_effect.length).times {|level1|
    mse_s[level1] = []
    v = pseudo_effect[level1]
#print "v = #{v} \n"
    (v.length).times {|factor1|
      mse_tmp = Array.new(v[factor1].length, 0)
      (v[factor1].length).times {|i|
        if (! v[check][i].nil? ) then 
          mse_tmp[i] = v[factor1][i] - v[check][i]
        else
          mse_tmp[i] = v[factor1][i]
        end
#print "i = #{i}, factor1 = #{factor1}, check = #{check}\n"
#print "v[factor1][i], v[check][i] = #{v[factor1][i]}, #{v[check][i]} \n"
#print "mse_tmp[i] = #{mse_tmp[i]}\n"
      }
      mse_s[level1][factor1] = mse_tmp
    }
  }
print "mse_s = #{mse_s}\n"

  mse = Array.new(mse_s[0].length)
  (mse_s[0].length).times {|i|
    q1 = Array.new(mse_s.length)
    (mse_s.length).times {|j|
      q1[j] = mse_s[j][i].flatten
      print "       i,j,q1[j] = #{i}, #{j}, #{q1[j]}\n"
    }
    str = 0
    k = 0
    0.upto(q1.length-2) { |i1|
#     must start at "1"
      1.upto(q1.length-i1-1) { |i2|
print "q1.len, i1, i2 = #{q1.length}, #{i1}, #{i2}\n"
        k1 = i1
        k2 = i2 + k1
        print "k1,k2 = #{k1}, #{k2}\n"
        wk = Array.new(q1[i1].length)
        0.upto(q1[i1].length-1) { |l|
          wk[l] = (q1[i1][l] - q1[i1+i2][l])**2
        }
        str += wk.average
        k += 1
print "str, k = #{str}, #{k}\n"
      }
      mse[i] = str/(k + 0.0)
print "mse[#{i}] = #{mse[i]}\n"
    }
  }
  mse
end
    
def calc_edge_information(mse1, threshold=1, zero_flg=nil, single_flg=nil)
  in1 = []
  out1 = []
  weight = []
  (mse1.length).times {|i|
    q = mse1[i]
    (q.length).times {|j|
 print "mse #{i} #{j}\n"
      if(mse1[i][j]>threshold) then
          in1 << i
          out1 << j
          weight << mse1[i][j]
      end
      if(i == j) then
      end
    }
  }
#p in1
#p out1
#p weight
  in1b = []
  out1b = []
  weightb = []
  (in1.length).times {|i|
    q = (out1.length).times.select{|j| in1[i] == out1[j]}
#print "q = #{q}\n"
    r = (q.length).times.select{|k| in1[q[k]] == out1[i]}
#print "r = #{r}\n"
    if ( r.length > 0) then
      in1b << in1[i]
      out1b << out1[i]
#      weightb = weight[i]
      weightb << weight[i]
    end
  }
  return [in1b, out1b, weightb]
end

def calc_node_size(data_set)
  node_size = Array.new(data_set.length)
  (data_set.length).times {|i|
    data_set_i = data_set[i][0]
    q = []
    (data_set_i.length).times {|j|
      q << data_set_i[j].flatten.average
    }
    node_size[i] = q.standard_devitation
  }
  node_size
end

def calc_node_color(results_confound)
  node_color = Array.new(results_confound.length)

  (results_confound.length).times {|i|
    result_i = results_confound[i]

#print "result_i[1] #{result_i[1]}\n"
    single = result_i[1].flatten.include?(TRUE)
#print "result_i[2] #{result_i[2]}\n"
    confound = result_i[2].flatten.include?(TRUE)
#print "single, confound = #{single} #{confound}\n"
    col_no = 1

    if (single) then col_no = 2 end
    if (confound) then col_no = 3 end

    node_color[i] = col_no
  }

  node_color
end

def calc_node_shape(results_confound)
  node_shape = Array.new(results_confound.length)
  (results_confound.length).times {|i|
    result_i = results_confound[i]
    q = []
    noise = result_i[0].flatten.include?(TRUE)
    shape_no = 1
    if (noise) then shape_no = 2 end
    node_shape[i] = shape_no
  }
  node_shape
end

def calc_effects_c_using_table(data_set_i_j, table_factor_check_is_i, mean_remove_flg=TRUE)
  table1 = table_factor_check_is_i

  factor_no = table1[0].length
  line_no = table1.length
  repeat_no = data_set_i_j.length
  
  if (mean_remove_flg == TRUE) then
    mean1 = (data_set_i_j.flatten).inject{|sum, n| sum + n} / repeat_no
  else
    mean1 = 0
  end

print "line_no = #{line_no}\n"

  ff = Array.new(factor_no, 0)
  factor_no.times {|i|
    nn = Array.new(line_no, 0)
    line_no.times {|p|
      nn[p] = table1[p][i]
    }

# variable "vvv" is not used.
print "nn.compact.max = #{nn.compact.max}\n"
p nn.compact.max

    a = Array.new((nn.compact.max).to_i, 0)
    b = Array.new((nn.compact.max).to_i, 0)

    line_no.times {|j|
      repeat_no.times {|k|
#p table1[j][i]
#p a
        a[table1[j][i].to_i-1] += data_set_i_j[k][j]
        b[table1[j][i].to_i-1] += 1.0
      }
    }

    if (mean_remove_flg == TRUE) then
print "a, b = #{a}, #{b}\n"
      wk = []
      (a.length).times {|m|
        if (b[m] != 0) then
          wk << a[m] / b[m]
        end
      }
      mean = wk.average
      wk2 = []
      (a.length).times {|m|
        if (b[m] != 0) then
          wk2 << a[m] / b[m] - mean
        end
      }
print "wk2 = #{wk2}\n"
      ff[i] = wk2
    else
      wk2 = []
      (a.length).times {|m|
        wk2 << a[m] / b[m]
      }
      ff[i] = wk2
    end
  }
print "ff = #{ff}\n"
  ff
end


def new_show_interaction_graph_all_c(mse1,data_set,results_confound,threshold=nil)
  if (threshold.nil?) then
print "******set threshold\n"
    threshold = 2 * total_sd(data_set)
  end

  edge_list_and_weight = calc_edge_information(mse1,threshold,zero_flg=nil, single_flg=nil)
#print "edge_list_and_weight = \n"
#p edge_list_and_weight

  node_size = calc_node_size(data_set)
#print "node_size\n"
#p node_size
  node_color = calc_node_color(results_confound)
#print "node_color\n"
#p node_color
  node_shape = calc_node_shape(results_confound)
#print "node_shape\n"
#p node_shape
  return [edge_list_and_weight, node_size, node_color, node_shape]
end

def total_sd(data_set)
  sd_mean = []
  (data_set.length).times {|i|
    (data_set[i][0].length).times {|j|
      (data_set[i][0][j][0].length).times {|k2|
        q = Array.new(data_set[i][0][j].length)
        (data_set[i][0][j].length).times {|k1|
          q[k1] = data_set[i][0][j][k1][k2]
        }
        sd_mean << q.standard_devitation
      }
    }
  }
  sd_mean.average
end

def new_arrange_the_confounding_all(results_confound)
  noise = []
  none = []
  single = []
  confound = []
#  print "results_confound.length = #{results_confound.length} \n"
#  0.upto(results_confound.length-1) {|i|
  (results_confound.length).times {|i|
    t0, t1, t2, t3 = new_arrange_the_confounding_i(results_confound,i)
    noise    << t0
    none     << t1
    single   << t2
    confound << t3
  }

#  open("result.out", "w") { |f|
    print "-----------Confounding Information Begin---------------\n"
    print "{\n"
    print "Noise:{ #{(noise - [""]).join(",")} },\n"
    print "None:{ #{(none - [""]).join(",")} },\n"
    print "Single:{ #{(single - [""]).join(",")} },\n"
    print "Confound:{ #{(confound - [""]).join(",")} }\n"
    print "}\n"
    print "-----------Confounding Information End---------------\n"
#  }
  return [noise,none,single,confound]
end

def new_arrange_the_confounding_i(results_confound,check)
  results_i = results_confound[check]
  noise = ""
  none = ""
  single = ""
  confound = ""
  flg = TRUE
  (results_i[0].length).times {|i|
    (results_i[0][0].length).times {|j|
      if (i>j) then
        if (results_i[0][i][j]) then
          if (flg) then
            noise = (check).to_s + ":["
            flg = FALSE
            noise += "[" + (i).to_s + "," + (j).to_s + "]"
          else
            noise += ",[" + (i).to_s + "," + (j).to_s + "]"
          end
        end
      end
    }
  }
  if (!flg) then
    noise[noise.length-1] = noise[noise.length-1].to_s + "]"
  end

  flg = TRUE
  (results_i[0].length).times {|i|
    (results_i[0][0].length).times {|j|
      if (i>j) then
        if (!results_i[1][i][j] & !results_i[2][i][j]) then
          if (flg) then
            none = (check).to_s + ":["
            flg = FALSE
            none += "[" + (i).to_s + "," + (j).to_s + "]"
          else
            none += ",[" + (i).to_s + "," + (j).to_s + "]"
          end
        end
      end
    }
  }
#print "none = #{none}\n"
  if (!flg) then
    none[none.length-1] = none[none.length-1].to_s + "]"
  end

  flg = TRUE
  (results_i[0].length).times {|i|
    (results_i[0][0].length).times {|j|
      if (i>j) then
        if (results_i[1][i][j] & !results_i[2][i][j]) then
          if (flg) then
            single = (check).to_s + ":["
            flg = FALSE
            single += "[" + (i).to_s + "," + (j).to_s + "]"
          else
            single += ",[" + (i).to_s + "," + (j).to_s + "]"
          end
        end
      end
    }
  }
  if (!flg) then
    single[single.length-1] = single[single.length-1].to_s + "]"
  end

  flg = TRUE
  (results_i[0].length).times {|i|
    (results_i[0][0].length).times {|j|
      if (i>j) then
        if (results_i[2][i][j]) then
          if (flg) then
            confound = (check).to_s + ":["
            flg = FALSE
            confound += "[" + (i).to_s + "," + (j).to_s + "]"
          else
            confound += ",[" + (i).to_s + "," + (j).to_s + "]"
          end
        end
      end
    }
  }
  if (!flg) then
    confound[confound.length-1] = confound[confound.length-1].to_s + "]"
  end

  return [noise, none, single, confound]
end

def arrange_the_network(network_information)
#  p network_information
#  open("result.out", "a") { |f|
    print "-----------Network Information Begin---------------\n"
    edge = network_information[0]
#  print "edge = #{edge}\n"
    print "{\n"
    print "Edge_Out_In_Weight:[\n"
    (edge[0].length).times {|i|
      if (i < edge[0].length-1) then
        print "[#{edge[1][i]}, #{edge[0][i]}, #{edge[2][i]}],\n"
      else
        print "[#{edge[1][i]}, #{edge[0][i]}, #{edge[2][i]}]],\n"
      end
    }
	
    size = network_information[1]
    print "Node_Size:["
    (size.length).times {|i|
      print ", " unless i == 0
      printf "%.4e",size[i]
    }
    print "],\n"

    color = network_information[2]
    print "Node_Color:["
    (color.length).times {|i|
      if (i != 0) then print "," end
      print "#{color[i]}"
    }
    print "],\n"

    shape = network_information[3]
    print "Node_Shape:["
    (shape.length).times {|i|
      if (i != 0) then print "," end
      print "#{shape[i]}"
    }
    print "]\n"
    print "}\n"
    print "-----------Network Information End---------------\n"
#  }
end


#### Start.
#open("_level_value.json") do |io|
#open("_level_value0.json") do |io|
#open("_level_value1.json") do |io|
#open("_level_value2.json") do |io|
#open("_level_value3.json") do |io|
#open("_level_value3-2.json") do |io|
#open("_level_value4.json") do |io|
open("_level_value4-2.json") do |io|
#open("_level_value4-3.json") do |io|
   @level = JSON.load(io)
end

lv = @level.values
print "lv = #{lv}\n"

vv2 = Marshal.load(Marshal.dump(lv))
# vv1 = Array.new(vv2.length)
# 0.upto(vv2.length-1) {|i| vv1[i] = [*1..vv2[i].length]}

nfactors = vv2.length
nlevels = Array.new(nfactors)
nfactors.times{|i| nlevels[i] = vv2[i].length}

vv1 = Array.new(vv2.length)
(vv2.length).times {|i| vv1[i] = [*1..vv2[i].length] }
p "vv1 = #{vv1}"
p "nfactors = #{nfactors}"
p "nlevels = #{nlevels}"

lvlstr = nlevels.join(",")
rstr ="library(\"DoE.base\")
  nlevels<-c(#{lvlstr})
  plan<-as.matrix(oa.design(nfactors=#{nfactors},nlevels=nlevels,seed=1))"

# print "rstr = #{rstr}\n"

r.eval_R rstr 

plan = r.plan
(plan.length).times {|i| print "#{plan[i]} \n" }

#fout_orth(plan)

repeat1 = 10
fun = method(:f)

### OACISの計算順をシミュレートする場合は、以下を有効にしてください。
#   (その場合は、Rと同等の処理をさせる部分をコメントにしてください。)

rev_orth, param_set = or_get_parameter_set(vv2, nfactors, plan)
#dump_array(rev_orth, nfactors, "rev_orth")
#dump_array(param_set, nfactors, "param_set")
#fout_revorth(rev_orth, nfactors)

nparam = count_param(param_set, nfactors)
print "nparam = #{nparam}\n"

sdata = or_sampling_data(fun, param_set, repeat1)
data_set = arrange_data_format(nfactors, nlevels, rev_orth, sdata, repeat1)

### OACISのシミュレーションを有効にするのはここまで。


### Rプログラムと同じ動作をさせる場合は以下のコメントを外してください。
#   (その場合は、OACISのシミュレート部分をコメントにしてください。)
=begin
data_set = sampling_data(fun, v1=vv1, v2=vv2, repeat1=repeat1, plan=plan)
=end
### Rと同等の動作をさせる場合にコメント解除するのはここまで。


#print "data_set = \n"
#dump_data_set(data_set, nfactors, nlevels, repeat1)
#fout_data_set(data_set, nfactors, nlevels, repeat1)

results_confound = finding_interaction_all_c(data_set)

print "results_confound =\n"
(results_confound.length).times {|i|
  (results_confound[i].length).times {|j|
    (results_confound[i][j].length).times {|k|
       print "#{results_confound[i][j][k]}\n"
    }
    print "\n"
  }
}

mse1 = check_interaction_all_c(data_set)

v = new_show_interaction_graph_all_c(mse1, data_set, results_confound, threshold=nil)

tmp1 = new_arrange_the_confounding_all(results_confound)

arrange_the_network(v)
