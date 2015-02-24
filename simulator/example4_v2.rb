#!/home/osaka/.rvm/rubies/ruby-1.9.3-p551/bin/ruby

require 'json'

def norm_rand(mu=0,sigma=1.0)
  r1,r2=Random.rand(),Random.rand();
  Math::sqrt(-1*Math::log(r1))*Math::cos(2*Math::PI*r2)*sigma+mu
end

# x0=ARGV[0].to_f
# x1=ARGV[1].to_f
# x2=ARGV[2].to_f
# x3=ARGV[3].to_f
# x4=ARGV[4].to_f
# x5=ARGV[5].to_f
# x6=ARGV[6].to_f

parsed = JSON.load(open('./_input.json'))

x0 = parsed["x0"]
x1 = parsed["x1"]
x2 = parsed["x2"]
x3 = parsed["x3"]
x4 = parsed["x4"]
x5 = parsed["x5"]
x6 = parsed["x6"]
x7 = parsed["x7"]
x8 = parsed["x8"]
x9 = parsed["x9"]

x10 = parsed["x10"]
x11 = parsed["x11"]
x12 = parsed["x12"]
x13 = parsed["x13"]
x14 = parsed["x14"]
x15 = parsed["x15"]
x16 = parsed["x16"]
x17 = parsed["x17"]
x18 = parsed["x18"]
x19 = parsed["x19"]

x20 = parsed["x20"]
x21 = parsed["x21"]
x22 = parsed["x22"]
x23 = parsed["x23"]
x24 = parsed["x24"]
x25 = parsed["x25"]
x26 = parsed["x26"]
x27 = parsed["x27"]
x28 = parsed["x28"]
x29 = parsed["x29"]

x30 = parsed["x30"]
x31 = parsed["x31"]
x32 = parsed["x32"]
x33 = parsed["x33"]
x34 = parsed["x34"]
x35 = parsed["x35"]
x36 = parsed["x36"]
x37 = parsed["x37"]
x38 = parsed["x38"]
x39 = parsed["x39"]

x40 = parsed["x40"]
x41 = parsed["x41"]
x42 = parsed["x42"]
x43 = parsed["x43"]
x44 = parsed["x44"]
x45 = parsed["x45"]
x46 = parsed["x46"]
x47 = parsed["x47"]
x48 = parsed["x48"]
x49 = parsed["x49"]

x50 = parsed["x50"]
x51 = parsed["x51"]
x52 = parsed["x52"]
x53 = parsed["x53"]
x54 = parsed["x54"]
x55 = parsed["x55"]
x56 = parsed["x56"]
x57 = parsed["x57"]
x58 = parsed["x58"]
x59 = parsed["x59"]

x60 = parsed["x60"]
x61 = parsed["x61"]
x62 = parsed["x62"]

     # Example 4 ver2

result = { "y" => 100*x30*100*x50+100*x26*100*x9+100*x4*100*x50+100*x35*100*x19+100*x62*100*x56+100*x53*100*x51+100*x29*100*x60+100*x58*100*x60+100*x7*100*x49+100*x10*100*x49+100*x49*100*x18+100*x40*100*x27+100*x0*100*x35+100*x47*100*x1+100*x34*100*x46+norm_rand(mu=0,sigma=1.0) }

JSON.dump(result, open('./_output.json','w'))
