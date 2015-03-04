# -*- coding: utf-8 -*-
require 'json'
require_relative '../../lib/OACIS_module.rb'
require_relative '../../lib/OACIS_module_data.rb'
require 'rsruby'


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


class IadModule < OacisModule

  def self.definition
    h = {}
    h["target_field"] = "y"
    h
  end


  def initialize(input_data)
    # get orthogonal table
    r = RSRuby.instance

    open("_level_value.json") do |io|
      @level = JSON.load(io)
    end

    lv = @level.values
#    print "lv = #{lv}\n"
    
    @vv2 = Marshal.load(Marshal.dump(lv))
    @nfactors = @vv2.length
    @nlevels = Array.new(@nfactors)
    @nfactors.times{|i| @nlevels[i] = @vv2[i].length}

    vv1 = Array.new(@vv2.length)
    (@vv2.length).times {|i| vv1[i] = [*1..@vv2[i].length] }
    p "vv1 = #{vv1}"
    p "nfactors = #{@nfactors}"
    p "nlevels = #{@nlevels}"

    lvlstr = @nlevels.join(",")
    rstr ="library(\"DoE.base\")
nlevels<-c(#{lvlstr})
plan<-as.matrix(oa.design(nfactors=#{@nfactors},nlevels=nlevels,seed=1))"

#   print "rstr = #{rstr}\n"
    r.eval_R rstr 
    @plan = r.plan
    fout_orth(@plan)

#    (plan.length).times {|i| print "#{plan[i]} \n" }
    @repeat1 = input_data["_target"]["RunsCount"]

    # create parameter set
    @rev_orth, @parameter_set = or_get_parameter_set(@vv2, @nfactors, @plan)

    # count number of parameter_sets.
    @nparam = count_param(@parameter_set, @nfactors)

    print "number of parameter sets = #{@nparam}\n"

    # create result array.
    @result_array = Array.new(@nparam, 0)

# binding.pry    
    super(input_data)
    @fin_flag = false
    @n_finished = 0
    @scope = 2000
  end


  private
  #override
  def generate_runs

    print "num_iterations = #{@num_iterations}\n"
    counter = 0
    @nregisterd = 0
    @nfactors.times {|k|
      (@parameter_set[k].length).times {|j|
#        print "param_set[#{k}][#{j}]\n"
        (@parameter_set[k][j].length).times {|i|
#          print "#{parameter_set[k][j][i]}\n"
          if (counter >= @n_finished and @n_finished + @scope > counter) then
            ret = module_data.set_input(@num_iterations, @nregisterd, @parameter_set[k][j][i])
            @nregisterd += 1
          end
          counter += 1
        }
      }
    }

#    @nregisterd = counter
    print "number of registerd parameter set = #{@nregisterd}\n"

# binding.pry
    super
  end


  #override
  def evaluate_runs
    super
    # get result of _output.json
    # result_array = Array.new(@nparam, 0)
    #(@nparam).times do |i| 
    (@nregisterd).times do |i| 
      y_array = module_data.get_output(@num_iterations, i)
      @result_array[i + @n_finished] = y_array
    end
    @n_finished += @nregisterd
    print "n_finished = #{@n_finished}\n"
    # get imulatior ID
    # sim = module_data.data["_input_data"]["_target"]["Simulator"]

    if @n_finished == @nparam then
      # main operation
      data_set = arrange_data_format(@nfactors, @nlevels, @rev_orth, @result_array, @repeat1)
      #debug
      #dump_data_set(data_set, @nfactors, @nlevels, @repeat1)
      fout_data_set(data_set, @nfactors, @nlevels, @repeat1)

      results_confound = finding_interaction_all_c(data_set)
      mse1 = check_interaction_all_c(data_set)
      v = new_show_interaction_graph_all_c(mse1, data_set, results_confound, threshold=nil)
      tmp1 = new_arrange_the_confounding_all(results_confound)
      arrange_the_network(v)
      @fin_flag=true
    end
  end


  #override
  def finished?
#binding.pry
    @fin_flag
  end

  #override
  def get_target_fields(result)
    result.try(:fetch, module_data.data["_input_data"]["target_field"])
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


  def or_make_data_set(rev_orth1, v2)
    param = Marshal.load(Marshal.dump(rev_orth1))
    #   print "param"
    #   p param
    #   print "v2"
    #   p v2
    (param.length).times {|k|
      (param[k].length).times {|j|
        (param[k][j].length).times {|i|
#          print "v2[ #{k} ] = #{v2[k]}\n"
#          print "param[#{k}][#{j}][#{i}] = #{param[k][j][i]}\n"
          param[k][j][i] = v2[i][param[k][j][i].to_i - 1]
        }
      }
    }
    #   p param
    param
  end


  def or_get_parameter_set(v2, nfactors=1, plan=nil)
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
#    dump_array(param_set, nfactors, "param_set")
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


  def arrange_data_format(nfactors, nlevels, rev_orth, data0, repeat)
    data_set = Array.new(nfactors).map{Array.new(2, 0)}
    num_experiment = rev_orth[0][0].length
    print "nfactors = #{nfactors}\n"
    print "nlevels = #{nlevels}\n"
    print "repeat = #{repeat}\n"
    print "num_experiment = #{num_experiment}\n"

    data1 = conv_flat_data(data0, num_experiment, nlevels)

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
        f.print "factor = #{factor}\n"
        (nlevels[factor]).times {|level|
          f.print "  level = #{level}\n"
          repeat.times {|nrep|
            f.print "    repeat #{nrep}  "
            f.print "#{data_set[factor][0][level][nrep]}\n"
          }
          f.print "    Rev orth\n"
          (data_set[factor][1][level].length).times {|nexp|
            f.print "      #{data_set[factor][1][level][nexp]}\n"
          }
        }
      }
    }
  end


  def fout_orth(plan)
    open("orth.out", "w") { |f|
      (plan.length).times {|i| f.print "#{plan[i]} \n" }
    }
  end


  def fout_revorth(rev_orth, nfactors)
    open("revorth.out", "w") { |f|
      nfactors.times {|k|
        (rev_orth[k].length).times {|j|
          f.print "[#{k}][#{j}]\n"
          (rev_orth[k][j].length).times {|i|
            f.print "#{rev_orth[k][j][i]}\n"
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

        print "sd1, sd2 = #{sd1}, #{sd2}\n"

        threshold = 2
        if [sd1, sd2].min != 0 then
          if ( (sd1 - sd2).abs/[sd1, sd2].min > threshold) then
            t[i][j] = TRUE
          else
            t[i][j] = FALSE
          end
        else
          if sd1 == sd2 then
            t[i][j] = FALSE
          else
            t[i][j] = TRUE
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
          print "answer TRUE: #{i}, #{j}\n"
        else
          t1[i][j] = FALSE
          print "answer FALSE: #{i}, #{j}\n"
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

#    print "q1 start\n"
#    p q1
#    p q2
#    print "q1 end\n"
    (q1[0].length).times {|kk|
      v1 = Array.new(q1.length, 0)
      v2 = Array.new(q1.length, 0)
      (q1.length).times {|k|
#        v1[k] = q1[k][kk].flatten
#        v2[k] = q2[k][kk].flatten
        v1[k] = q1[k][kk]
        v2[k] = q2[k][kk]
      }
#      print "v1 #{v1}\n"
#      print "v2 #{v2}\n"
      vv1[kk] = v1.average
      vv2[kk] = v2.average
    }

    v3 = Array.new(vv1.length, 0)
    (vv1.length).times {|i|
      v3[i] = (vv1[i] - vv2[i]).abs
    }
#    print "vv1 #{vv1}\n"
#    print "vv2 #{vv2}\n"
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
#      table_factor_check_is_i = as_list_orthogonal_table(table_list[level].transpose)
#      print "table_list[level] = #{table_list[level]}\n"
#      table_factor_check_is_i = table_list[level].transpose
      table_factor_check_is_i = table_list[level]
#      print "table_factor_check_is_i = #{table_factor_check_is_i}\n"
      data_set_i_j = data_set_i[0][level]
      v = calc_effects_c_using_table(data_set_i_j, table_factor_check_is_i)
#      print "v = #{v}\n"
      pseudo_effect[level] = v
    }

    print "pseudo_effect = #{pseudo_effect}\n"

    mse_s = Array.new(pseudo_effect.length, 0)

    (pseudo_effect.length).times {|level1|
      mse_s[level1] = []
      v = pseudo_effect[level1]
#      print "v = #{v} \n"
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
#    print "mse_s = #{mse_s}\n"

    mse = Array.new(mse_s[0].length)
    (mse_s[0].length).times {|i|
      q1 = Array.new(mse_s.length)
      (mse_s.length).times {|j|
        q1[j] = mse_s[j][i].flatten
#        print "       i,j,q1[j] = #{i}, #{j}, #{q1[j]}\n"
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
#      #print "q = #{q}\n"
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

#    print "line_no = #{line_no}\n"

    ff = Array.new(factor_no, 0)
    factor_no.times {|i|
      nn = Array.new(line_no, 0)
      line_no.times {|p|
        nn[p] = table1[p][i]
      }

      # variable "vvv" is not used.
#      print "nn.compact.max = #{nn.compact.max}\n"
#      p nn.compact.max

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
#        print "a, b = #{a}, #{b}\n"
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
#        print "wk2 = #{wk2}\n"
        ff[i] = wk2
      else
        wk2 = []
        (a.length).times {|m|
          wk2 << a[m] / b[m]
        }
        ff[i] = wk2
      end
    }
#    print "ff = #{ff}\n"
    ff
  end


  def new_show_interaction_graph_all_c(mse1,data_set,results_confound,threshold=nil)
    if (threshold.nil?) then
      print "******set threshold\n"
      threshold = 2 * total_sd(data_set)
    end

    edge_list_and_weight = calc_edge_information(mse1,threshold,zero_flg=nil, single_flg=nil)
    print "edge_list_and_weight = \n"
    p edge_list_and_weight

    node_size = calc_node_size(data_set)
    print "node_size\n"
    p node_size
    node_color = calc_node_color(results_confound)
    print "node_color\n"
    p node_color
    node_shape = calc_node_shape(results_confound)
    print "node_shape\n"
    p node_shape
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

    print "results_confound.length = #{results_confound.length} \n"
    #  0.upto(results_confound.length-1) {|i|
    (results_confound.length).times {|i|
      t0, t1, t2, t3 = new_arrange_the_confounding_i(results_confound,i)
      noise    << t0
      none     << t1
      single   << t2
      confound << t3
    }

    resultstr =  "-----------Confounding Information Begin---------------\n"
    resultstr += "{\n"
    resultstr += "Noise:{ #{(noise - [""]).join(",")} },\n"
    resultstr += "None:{ #{(none - [""]).join(",")} },\n"
    resultstr += "Single:{ #{(single - [""]).join(",")} },\n"
    resultstr += "Confound:{ #{(confound - [""]).join(",")} }\n"
    resultstr += "}\n"
    resultstr += "-----------Confounding Information End---------------\n"

    # print result.
    print resultstr

    # save result to file.
    open("result.out", "w") { |f|
      f.print resultstr
    }
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
    print "none = #{none}\n"
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
    resultstr =  "-----------Network Information Begin---------------\n"
    edge = network_information[0]
    #  print "edge = #{edge}\n"
    resultstr += "{\n"
    resultstr += "Edge_Out_In_Weight:[\n"
    (edge[0].length).times {|i|
      if (i < edge[0].length-1) then
        resultstr += "[#{edge[1][i]}, #{edge[0][i]}, #{edge[2][i]}],\n"
      else
        resultstr += "[#{edge[1][i]}, #{edge[0][i]}, #{edge[2][i]}]],\n"
      end
    }
	
    size = network_information[1]
    resultstr += "Node_Size:["
    (size.length).times {|i|
      resultstr += ", " unless i == 0
      resultstr += sprintf "%.4e",size[i]
    }
    resultstr += "],\n"

    color = network_information[2]
    resultstr += "Node_Color:["
    (color.length).times {|i|
      if (i != 0) then resultstr += "," end
      resultstr += "#{color[i]}"
    }
    resultstr += "],\n"

    shape = network_information[3]
    resultstr += "Node_Shape:["
    (shape.length).times {|i|
      if (i != 0) then resultstr += "," end
      resultstr += "#{shape[i]}"
    }
    resultstr += "]\n"
    resultstr += "}\n"
    resultstr += "-----------Network Information End---------------\n"

    # print result.
    print resultstr

    # save result to file.
    open("result.out", "a") { |f|
      f.print resultstr
    }

  end

  ### End of function definitions.
end
