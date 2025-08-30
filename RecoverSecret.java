import com.google.gson.*;
import com.google.gson.JsonParser;
import java.io.*;
import java.math.BigInteger;
import java.util.*;

/**
 * RecoverSecre
 * Reads JSON from input.json. Finds secret (f(0)) using Lagrange interpolation on k sized subsets,
 * picks the polynomial that matches the most shares, prints the secret and the wrong share keys.
 *
 * Output format:
 * Secret: <decimal-string>
 * Wrong shares: <comma-separated list of keys>   (or "None" if all match)
 */
public class RecoverSecret {
    // Simple rational number class using BigInteger
    static final class Rational {
        BigInteger num;
        BigInteger den; // always > 0

        Rational(BigInteger n, BigInteger d) {
            if (d.signum() == 0) throw new ArithmeticException("zero denominator");
            if (d.signum() < 0) { n = n.negate(); d = d.negate(); }
            BigInteger g = n.gcd(d);
            if (!g.equals(BigInteger.ONE)) {
                n = n.divide(g);
                d = d.divide(g);
            }
            this.num = n;
            this.den = d;
        }
        Rational(long a) { this(BigInteger.valueOf(a), BigInteger.ONE); }
        Rational(BigInteger a) { this(a, BigInteger.ONE); }

        static Rational zero() { return new Rational(BigInteger.ZERO, BigInteger.ONE); }
        static Rational one() { return new Rational(BigInteger.ONE, BigInteger.ONE); }

        Rational add(Rational o) {
            BigInteger n = this.num.multiply(o.den).add(o.num.multiply(this.den));
            BigInteger d = this.den.multiply(o.den);
            return new Rational(n, d);
        }
        Rational sub(Rational o) {
            BigInteger n = this.num.multiply(o.den).subtract(o.num.multiply(this.den));
            BigInteger d = this.den.multiply(o.den);
            return new Rational(n, d);
        }
        Rational mul(Rational o) {
            return new Rational(this.num.multiply(o.num), this.den.multiply(o.den));
        }
        Rational div(Rational o) {
            if (o.num.signum() == 0) throw new ArithmeticException("div by zero");
            return new Rational(this.num.multiply(o.den), this.den.multiply(o.num));
        }
        BigInteger toBigIntegerIfIntegral() {
            if (den.equals(BigInteger.ONE)) return num;
            if (num.mod(den).equals(BigInteger.ZERO)) return num.divide(den);
            return null;
        }
        @Override public String toString() {
            if (den.equals(BigInteger.ONE)) return num.toString();
            return num.toString() + "/" + den.toString();
        }
    }

    static final class Share {
        BigInteger x;
        BigInteger y;
        String key;
        Share(String key, BigInteger x, BigInteger y) { this.key = key; this.x = x; this.y = y; }
    }

    static Rational lagrangeF0(List<Share> subset) {
        int k = subset.size();
        Rational acc = Rational.zero();
        for (int i = 0; i < k; ++i) {
            Share si = subset.get(i);
            Rational term = new Rational(si.y);
            Rational li_at0 = Rational.one();
            for (int j = 0; j < k; ++j) {
                if (j == i) continue;
                Share sj = subset.get(j);
                Rational numerator = new Rational(sj.x.negate());
                Rational denominator = new Rational(si.x.subtract(sj.x));
                li_at0 = li_at0.mul(numerator.div(denominator));
            }
            term = term.mul(li_at0);
            acc = acc.add(term);
        }
        return acc;
    }

    static Rational lagrangeEvaluateAt(List<Share> subset, BigInteger x) {
        int k = subset.size();
        Rational acc = Rational.zero();
        for (int i = 0; i < k; ++i) {
            Share si = subset.get(i);
            Rational term = new Rational(si.y);
            Rational li = Rational.one();
            for (int j = 0; j < k; ++j) {
                if (j == i) continue;
                Share sj = subset.get(j);
                Rational numerator = new Rational(x.subtract(sj.x));
                Rational denominator = new Rational(si.x.subtract(sj.x));
                li = li.mul(numerator.div(denominator));
            }
            term = term.mul(li);
            acc = acc.add(term);
        }
        return acc;
    }

    static void combinations(int n, int k, int start, ArrayList<Integer> cur, List<List<Integer>> out) {
        if (cur.size() == k) { out.add(new ArrayList<>(cur)); return; }
        for (int i = start; i <= n - (k - cur.size()); ++i) {
            cur.add(i);
            combinations(n, k, i + 1, cur, out);
            cur.remove(cur.size() - 1);
        }
    }

    public static void main(String[] args) throws Exception {
        // 🔹 Instead of stdin, read from input.json
        String jsonText = new String(java.nio.file.Files.readAllBytes(java.nio.file.Paths.get("input.json")));
        jsonText = jsonText.trim();
        if (jsonText.isEmpty()) {
            System.err.println("No JSON input found in input.json");
            return;
        }

        Gson gson = new Gson();
        JsonObject root = JsonParser.parseString(jsonText).getAsJsonObject();

        JsonObject keys = root.getAsJsonObject("keys");
        int n = keys.get("n").getAsInt();
        int k = keys.get("k").getAsInt();

        List<Share> shares = new ArrayList<>();
        for (Map.Entry<String, JsonElement> entry : root.entrySet()) {
            String key = entry.getKey();
            if (key.equals("keys")) continue;
            JsonObject shareObj = entry.getValue().getAsJsonObject();
            int base = Integer.parseInt(shareObj.get("base").getAsString());
            String valueStr = shareObj.get("value").getAsString();
            BigInteger y;
            try { y = new BigInteger(valueStr, base); }
            catch (Exception ex) { y = new BigInteger(valueStr.toLowerCase(), base); }
            BigInteger x = new BigInteger(key);
            shares.add(new Share(key, x, y));
        }

        if (shares.size() != n) {
            System.err.println("Warning: actual shares count (" + shares.size() + ") != n (" + n + ")");
            n = shares.size();
        }

        List<List<Integer>> combs = new ArrayList<>();
        combinations(n, k, 0, new ArrayList<>(), combs);

        int bestMatches = -1;
        Rational bestSecret = null;
        List<Boolean> bestMask = null;

        for (List<Integer> comb : combs) {
            List<Share> subset = new ArrayList<>();
            for (int idx : comb) subset.add(shares.get(idx));
            Rational f0;
            try { f0 = lagrangeF0(subset); }
            catch (Exception ex) { continue; }

            int matches = 0;
            List<Boolean> mask = new ArrayList<>();
            for (Share s : shares) {
                Rational eval = lagrangeEvaluateAt(subset, s.x);
                BigInteger maybeInt = eval.toBigIntegerIfIntegral();
                boolean equal = (maybeInt != null && maybeInt.equals(s.y));
                mask.add(equal);
                if (equal) matches++;
            }
            if (matches > bestMatches) {
                bestMatches = matches;
                bestSecret = f0;
                bestMask = mask;
            }
        }

        if (bestSecret == null) {
            System.err.println("No valid polynomial found.");
            return;
        }

        BigInteger secretInt = bestSecret.toBigIntegerIfIntegral();
        System.out.println("Secret: " + (secretInt != null ? secretInt.toString() : bestSecret.toString()));

        List<String> wrongKeys = new ArrayList<>();
        for (int i = 0; i < shares.size(); ++i) {
            if (!bestMask.get(i)) wrongKeys.add(shares.get(i).key);
        }
        if (wrongKeys.isEmpty()) {
            System.out.println("Wrong shares: None");
        } else {
            System.out.println("Wrong shares: " + String.join(",", wrongKeys));
        }
    }
}



