package trd;
import java.util.Random;
public class Main {

    public static void main(String[] args) {
        for (int i = 0; i < 5; i++) {
            afficherFruit();
        }
    }

    public static void afficherFruit() {

        Random random = new Random();
        int nombre = random.nextInt(10); // 0 à 9

        String[] name = {
            "Pomme", "Banane", "Orange", "Fraise", "Mangue",
            "Kiwi", "Ananas", "Cerise", "Poire", "Raisin"
        };
//teddtt
        switch (nombre) {
            case 0:
                System.out.println("Hello " + name[nombre]);
                break;
            case 1:
                System.out.println("Hello " + name[nombre]);
                break;
            case 2:
                System.out.println("Hello " + name[nombre]);
                break;
            case 3:
                System.out.println("Hello " + name[nombre]);
                break;
            case 4:
                System.out.println("Hello " + name[nombre]);
                break;
            case 5:
                System.out.println("Hello " + name[nombre]);
                break;
            case 6:
                System.out.println("Hello " + name[nombre]);
                break;
            case 7:
                System.out.println("Hello " + name[nombre]);
                break;
            case 8:
                System.out.println("Hello " + name[nombre]);
                break;
            case 9:
                System.out.println("Hello " + name[nombre]);
                break;
            default:
                System.out.println("Valeur invalide");
        }
    }
}
